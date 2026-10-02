#include "placamera/dataset_formats.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <optional>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace placamera
{
    namespace
    {

        std::string trim(std::string text)
        {
            const auto first = text.find_first_not_of(" \t\n\f\v");
            if (first == std::string::npos)
            {
                return {};
            }
            const auto last = text.find_last_not_of(" \t\n\f\v");
            return text.substr(first, last - first + 1);
        }

        std::string lower(std::string text)
        {
            std::transform(text.begin(),
                           text.end(),
                           text.begin(),
                           [](unsigned char character) { return static_cast<char>(std::tolower(character)); });
            return text;
        }

        void appendUtf8(std::string& output, unsigned int codepoint)
        {
            if (codepoint <= 0x7F)
            {
                output.push_back(static_cast<char>(codepoint));
            }
            else if (codepoint <= 0x7FF)
            {
                output.push_back(static_cast<char>(0xC0 | (codepoint >> 6)));
                output.push_back(static_cast<char>(0x80 | (codepoint & 0x3F)));
            }
            else if (codepoint <= 0xFFFF)
            {
                output.push_back(static_cast<char>(0xE0 | (codepoint >> 12)));
                output.push_back(static_cast<char>(0x80 | ((codepoint >> 6) & 0x3F)));
                output.push_back(static_cast<char>(0x80 | (codepoint & 0x3F)));
            }
            else
            {
                output.push_back(static_cast<char>(0xF0 | (codepoint >> 18)));
                output.push_back(static_cast<char>(0x80 | ((codepoint >> 12) & 0x3F)));
                output.push_back(static_cast<char>(0x80 | ((codepoint >> 6) & 0x3F)));
                output.push_back(static_cast<char>(0x80 | (codepoint & 0x3F)));
            }
        }

        bool validUtf8(std::string_view text)
        {
            for (std::size_t index = 0; index < text.size();)
            {
                const auto lead = static_cast<unsigned char>(text[index]);
                if (lead < 0x80)
                {
                    ++index;
                    continue;
                }
                int length = 0;
                unsigned int codepoint = 0;
                if (lead >= 0xC2 && lead <= 0xDF)
                {
                    length = 2;
                    codepoint = lead & 0x1F;
                }
                else if (lead >= 0xE0 && lead <= 0xEF)
                {
                    length = 3;
                    codepoint = lead & 0x0F;
                }
                else if (lead >= 0xF0 && lead <= 0xF4)
                {
                    length = 4;
                    codepoint = lead & 0x07;
                }
                else
                {
                    return false;
                }
                if (index + static_cast<std::size_t>(length) > text.size())
                {
                    return false;
                }
                for (int offset = 1; offset < length; ++offset)
                {
                    const auto byte = static_cast<unsigned char>(text[index + static_cast<std::size_t>(offset)]);
                    if ((byte & 0xC0) != 0x80)
                    {
                        return false;
                    }
                    codepoint = (codepoint << 6) | (byte & 0x3F);
                }
                if ((length == 3 && (codepoint < 0x800 || (codepoint >= 0xD800 && codepoint <= 0xDFFF))) ||
                    (length == 4 && (codepoint < 0x10000 || codepoint > 0x10FFFF)))
                {
                    return false;
                }
                index += static_cast<std::size_t>(length);
            }
            return true;
        }

        std::string normalizeMetashapeXml(std::string_view input)
        {
            if (input.starts_with("\xEF\xBB\xBF"))
            {
                input.remove_prefix(3);
            }
            if (input.starts_with("\xFF\xFE") || input.starts_with("\xFE\xFF"))
            {
                const bool big_endian = input.starts_with("\xFE\xFF");
                input.remove_prefix(2);
                if (input.size() % 2 != 0)
                {
                    throw std::runtime_error("Metashape XML has an incomplete UTF-16 code unit");
                }
                std::string output;
                output.reserve(input.size() / 2);
                for (std::size_t index = 0; index < input.size(); index += 2)
                {
                    const auto first = static_cast<unsigned char>(input[index]);
                    const auto second = static_cast<unsigned char>(input[index + 1]);
                    const unsigned int unit = big_endian ? (first << 8) | second : (second << 8) | first;
                    unsigned int codepoint = unit;
                    if (unit >= 0xD800 && unit <= 0xDBFF)
                    {
                        if (index + 3 >= input.size())
                        {
                            throw std::runtime_error("Metashape XML has an unpaired UTF-16 high surrogate");
                        }
                        const auto next_first = static_cast<unsigned char>(input[index + 2]);
                        const auto next_second = static_cast<unsigned char>(input[index + 3]);
                        const unsigned int next =
                            big_endian ? (next_first << 8) | next_second : (next_second << 8) | next_first;
                        if (next < 0xDC00 || next > 0xDFFF)
                        {
                            throw std::runtime_error("Metashape XML has an invalid UTF-16 surrogate pair");
                        }
                        codepoint = 0x10000 + ((unit - 0xD800) << 10) + (next - 0xDC00);
                        index += 2;
                    }
                    else if (unit >= 0xDC00 && unit <= 0xDFFF)
                    {
                        throw std::runtime_error("Metashape XML has an unpaired UTF-16 low surrogate");
                    }
                    appendUtf8(output, codepoint);
                }
                return output;
            }
            if (!validUtf8(input))
            {
                throw std::runtime_error("Metashape XML must be UTF-8 or UTF-16 with a byte-order mark");
            }
            return std::string(input);
        }

        std::optional<std::string> xmlAttribute(const std::string& xml, const std::string& name);

        double finiteNumber(std::string_view text, std::string_view field)
        {
            const std::string normalized = trim(std::string(text));
            if (normalized.empty())
            {
                throw std::runtime_error("Metashape " + std::string(field) + " is empty");
            }
            std::size_t consumed = 0;
            double value = 0.0;
            try
            {
                value = std::stod(normalized, &consumed);
            }
            catch (const std::exception&)
            {
                throw std::runtime_error("Metashape " + std::string(field) + " is not a valid number: " + normalized);
            }
            if (consumed != normalized.size() || !std::isfinite(value))
            {
                throw std::runtime_error("Metashape " + std::string(field) + " must be a finite number: " + normalized);
            }
            return value;
        }

        int integer(std::string_view text, std::string_view field)
        {
            const std::string normalized = trim(std::string(text));
            if (normalized.empty())
            {
                throw std::runtime_error("Metashape " + std::string(field) + " is empty");
            }
            std::size_t consumed = 0;
            int value = 0;
            try
            {
                value = std::stoi(normalized, &consumed);
            }
            catch (const std::exception&)
            {
                throw std::runtime_error("Metashape " + std::string(field) + " is not a valid integer: " + normalized);
            }
            if (consumed != normalized.size())
            {
                throw std::runtime_error("Metashape " + std::string(field) + " must be an integer: " + normalized);
            }
            return value;
        }

        FrameProjectionModel metashapeProjection(const std::string& sensor)
        {
            const std::string value =
                lower(xmlAttribute(sensor, "type").value_or(xmlAttribute(sensor, "projection").value_or("frame")));
            if (value == "frame" || value == "perspective" || value == "pinhole")
            {
                return FrameProjectionModel::Perspective;
            }
            if (value == "fisheye")
            {
                return FrameProjectionModel::Fisheye;
            }
            if (value == "equidistant" || value == "equidistant_fisheye")
            {
                return FrameProjectionModel::EquidistantFisheye;
            }
            if (value == "equisolid" || value == "equisolid_fisheye")
            {
                return FrameProjectionModel::EquisolidFisheye;
            }
            if (value == "spherical")
            {
                return FrameProjectionModel::Spherical;
            }
            if (value == "cylindrical")
            {
                return FrameProjectionModel::Cylindrical;
            }
            throw std::runtime_error("Metashape sensor has unsupported projection type: " + value);
        }

        RollingShutterMode metashapeRollingShutter(const std::string& sensor)
        {
            RollingShutterMode mode = RollingShutterMode::Disabled;
            const std::regex property_pattern("<property\\b[^>]*/?>");
            for (auto it = std::sregex_iterator(sensor.begin(), sensor.end(), property_pattern);
                 it != std::sregex_iterator();
                 ++it)
            {
                const std::string property = it->str();
                const auto name = xmlAttribute(property, "name");
                const auto value = xmlAttribute(property, "value");
                if (!name || !value)
                {
                    continue;
                }
                if (*name == "rolling_shutter")
                {
                    const std::string enabled = lower(*value);
                    if (enabled == "true" || enabled == "1" || enabled == "yes")
                    {
                        mode = RollingShutterMode::Full;
                    }
                    else if (enabled != "false" && enabled != "0" && enabled != "no")
                    {
                        throw std::runtime_error("Metashape rolling_shutter property is invalid: " + *value);
                    }
                }
                else if (*name == "rolling_shutter_flags")
                {
                    mode = integer(*value, "rolling_shutter_flags") == 3 ? RollingShutterMode::Regularized
                                                                           : RollingShutterMode::Full;
                }
            }
            return mode;
        }

        std::vector<std::string> dataLines(std::istream& input)
        {
            std::vector<std::string> lines;
            std::string line;
            while (std::getline(input, line))
            {
                line = trim(std::move(line));
                if (!line.empty() && line.front() != '#')
                {
                    lines.push_back(std::move(line));
                }
            }
            return lines;
        }

        std::vector<double> numericLine(const std::string& line, std::size_t count)
        {
            std::istringstream input(line);
            std::vector<double> values;
            double value = 0.0;
            while (input >> value)
            {
                values.push_back(value);
            }
            if (values.size() != count)
            {
                throw std::runtime_error("numeric line has the wrong field count: " + line);
            }
            return values;
        }

        std::array<double, 9> transpose(const std::array<double, 9>& matrix)
        {
            return {matrix[0], matrix[3], matrix[6], matrix[1], matrix[4], matrix[7], matrix[2], matrix[5], matrix[8]};
        }

        ImportedCamera middleburyLine(const std::string& line)
        {
            std::istringstream input(line);
            ImportedCamera camera;
            if (!(input >> camera.imageName))
            {
                throw std::runtime_error("Middlebury camera line has no image name: " + line);
            }
            std::vector<double> values;
            double value = 0.0;
            while (input >> value)
            {
                values.push_back(value);
            }
            if (values.size() != 21)
            {
                throw std::runtime_error("Middlebury camera line requires 21 numeric fields: " + line);
            }
            std::copy_n(values.begin(), 9, camera.calibration.intrinsicMatrix.begin());
            std::array<double, 9> world_to_camera{};
            std::copy_n(values.begin() + 9, 9, world_to_camera.begin());
            camera.cameraToWorldRotation = transpose(world_to_camera);
            for (std::size_t row = 0; row < 3; ++row)
            {
                camera.center[row] = -(camera.cameraToWorldRotation[row * 3] * values[18] +
                                       camera.cameraToWorldRotation[row * 3 + 1] * values[19] +
                                       camera.cameraToWorldRotation[row * 3 + 2] * values[20]);
            }
            return camera;
        }

        std::vector<std::string> xmlBlocks(const std::string& xml, const std::string& tag)
        {
            std::vector<std::string> blocks;
            const std::regex pattern("<" + tag + "\\b[^>]*>[\\s\\S]*?</" + tag + ">");
            for (auto it = std::sregex_iterator(xml.begin(), xml.end(), pattern); it != std::sregex_iterator(); ++it)
            {
                blocks.push_back(it->str());
            }
            return blocks;
        }

        std::optional<std::string> xmlElementText(const std::string& xml, const std::string& tag)
        {
            const std::regex pattern("<" + tag + "\\b[^>]*>([\\s\\S]*?)</" + tag + ">");
            std::smatch match;
            if (!std::regex_search(xml, match, pattern) || match.size() < 2)
            {
                return std::nullopt;
            }
            return trim(match[1].str());
        }

        std::optional<std::string> xmlAttribute(const std::string& xml, const std::string& name)
        {
            const std::regex pattern("(?:^|[\\s<])" + name + "\\s*=\\s*[\"']([^\"']*)[\"']");
            std::smatch match;
            if (!std::regex_search(xml, match, pattern) || match.size() < 2)
            {
                return std::nullopt;
            }
            return match[1].str();
        }

        double xmlDoubleOr(const std::string& xml, const std::string& tag, double fallback)
        {
            const auto value = xmlElementText(xml, tag);
            return value.has_value() && !value->empty() ? finiteNumber(*value, "<" + tag + ">") : fallback;
        }

        struct MetashapeSensor
        {
            ImportedPixelCalibration calibration;
            RollingShutterMode rollingShutterMode = RollingShutterMode::Disabled;
        };

        std::unordered_map<int, MetashapeSensor> metashapeSensors(const std::string& xml)
        {
            std::unordered_map<int, MetashapeSensor> sensors;
            for (const std::string& sensor_block : xmlBlocks(xml, "sensor"))
            {
                const auto id = xmlAttribute(sensor_block, "id");
                if (!id.has_value())
                {
                    continue;
                }
                const std::regex resolution_pattern("<resolution\\b[^>]*/?>");
                std::smatch resolution;
                if (!std::regex_search(sensor_block, resolution, resolution_pattern))
                {
                    throw std::runtime_error("Metashape sensor has no resolution");
                }
                const auto width_text = xmlAttribute(resolution.str(), "width");
                const auto height_text = xmlAttribute(resolution.str(), "height");
                if (!width_text.has_value() || !height_text.has_value())
                {
                    throw std::runtime_error("Metashape resolution has no width or height");
                }
                const double width = finiteNumber(*width_text, "resolution width");
                const double height = finiteNumber(*height_text, "resolution height");
                if (!(width > 0.0) || !(height > 0.0))
                {
                    throw std::runtime_error("Metashape resolution width and height must be positive");
                }

                const auto calibrations = xmlBlocks(sensor_block, "calibration");
                if (calibrations.empty())
                {
                    throw std::runtime_error("Metashape sensor has no calibration");
                }
                std::string calibration = calibrations.front();
                for (const auto& block : calibrations)
                {
                    if (xmlAttribute(block, "class") == "adjusted")
                    {
                        calibration = block;
                        break;
                    }
                }
                const double focal = xmlDoubleOr(calibration, "f", 0.0);
                const double b1 = xmlDoubleOr(calibration, "b1", 0.0);
                const double b2 = xmlDoubleOr(calibration, "b2", 0.0);
                const double fx = xmlDoubleOr(calibration, "fx", focal + b1);
                const double fy = xmlDoubleOr(calibration, "fy", focal);
                if (fx <= 0.0 || fy <= 0.0)
                {
                    throw std::runtime_error("Metashape calibration has no valid focal length");
                }
                const double cx = xmlDoubleOr(calibration, "cx", 0.0);
                const double cy = xmlDoubleOr(calibration, "cy", 0.0);
                MetashapeSensor sensor;
                sensor.calibration.projectionModel = metashapeProjection(sensor_block);
                sensor.rollingShutterMode = metashapeRollingShutter(sensor_block);
                sensor.calibration.intrinsicMatrix = {
                    fx, b2, width * 0.5 + cx, 0.0, fy, height * 0.5 + cy, 0.0, 0.0, 1.0};
                auto& distortion = sensor.calibration.distortion;
                distortion.radialK1 = xmlDoubleOr(calibration, "k1", 0.0);
                distortion.radialK2 = xmlDoubleOr(calibration, "k2", 0.0);
                distortion.radialK3 = xmlDoubleOr(calibration, "k3", 0.0);
                distortion.radialK4 = xmlDoubleOr(calibration, "k4", 0.0);
                distortion.tangentialP1 = xmlDoubleOr(calibration, "p1", 0.0);
                distortion.tangentialP2 = xmlDoubleOr(calibration, "p2", 0.0);
                distortion.tangentialP3 = xmlDoubleOr(calibration, "p3", 0.0);
                distortion.tangentialP4 = xmlDoubleOr(calibration, "p4", 0.0);
                distortion.tangentialConvention = BrownTangentialConvention::Metashape;
                MetashapeCalibration exact;
                exact.f = fy;
                exact.cx = width * 0.5 + cx;
                exact.cy = height * 0.5 + cy;
                exact.b1 = fx - fy;
                exact.b2 = b2;
                exact.k1 = distortion.radialK1;
                exact.k2 = distortion.radialK2;
                exact.k3 = distortion.radialK3;
                exact.k4 = distortion.radialK4;
                exact.p1 = distortion.tangentialP1;
                exact.p2 = distortion.tangentialP2;
                exact.p3 = distortion.tangentialP3;
                exact.p4 = distortion.tangentialP4;
                exact.principalPointDecomposition = PrincipalPointDecomposition{width * 0.5, height * 0.5, cx, cy};
                sensor.calibration.metashapeCalibration = exact;
                sensors.emplace(integer(*id, "sensor id"), sensor);
            }
            if (sensors.empty())
            {
                throw std::runtime_error("Metashape document has no usable sensor");
            }
            return sensors;
        }

    } // namespace

    Result<std::vector<ImportedCamera>> readMiddleburyPar(std::istream& input)
    {
        try
        {
            auto lines = dataLines(input);
            if (!lines.empty())
            {
                std::istringstream first(lines.front());
                int count = 0;
                std::string trailing;
                if ((first >> count) && !(first >> trailing))
                {
                    lines.erase(lines.begin());
                }
            }
            if (lines.empty())
            {
                throw std::runtime_error("Middlebury par file has no camera records");
            }
            std::vector<ImportedCamera> cameras;
            cameras.reserve(lines.size());
            for (const auto& line : lines)
            {
                cameras.push_back(middleburyLine(line));
            }
            return Result<std::vector<ImportedCamera>>::success(std::move(cameras));
        }
        catch (const std::exception& error)
        {
            return Result<std::vector<ImportedCamera>>::failure(
                CameraErrorCode::ParseFailure, error.what(), "Middlebury par");
        }
    }

    Result<ImportedCamera> readEpflCamera(std::istream& input)
    {
        try
        {
            const auto lines = dataLines(input);
            if (lines.size() < 8)
            {
                throw std::runtime_error("EPFL .camera file requires at least eight numeric lines");
            }
            ImportedCamera camera;
            for (std::size_t row = 0; row < 3; ++row)
            {
                const auto values = numericLine(lines[row], 3);
                std::copy(values.begin(), values.end(), camera.calibration.intrinsicMatrix.begin() + row * 3);
            }
            const auto radial = numericLine(lines[3], 3);
            for (std::size_t index = 0; index < radial.size(); ++index)
            {
                camera.compatibility.sourceOnlyTerms.push_back(
                    {"epfl_radial_" + std::to_string(index + 1), radial[index]});
            }
            for (std::size_t row = 0; row < 3; ++row)
            {
                const auto values = numericLine(lines[row + 4], 3);
                std::copy(values.begin(), values.end(), camera.cameraToWorldRotation.begin() + row * 3);
            }
            const auto center = numericLine(lines[7], 3);
            std::copy(center.begin(), center.end(), camera.center.begin());
            return Result<ImportedCamera>::success(std::move(camera));
        }
        catch (const std::exception& error)
        {
            return Result<ImportedCamera>::failure(CameraErrorCode::ParseFailure, error.what(), "EPFL .camera");
        }
    }

    Result<std::vector<ImportedCamera>> parseMetashapeDocument(std::string_view xml)
    {
        try
        {
            const std::string document = normalizeMetashapeXml(xml);
            const auto sensors = metashapeSensors(document);
            std::vector<ImportedCamera> cameras;
            for (const std::string& block : xmlBlocks(document, "camera"))
            {
                const auto label = xmlAttribute(block, "label");
                const auto sensor_id = xmlAttribute(block, "sensor_id");
                const auto transform = xmlElementText(block, "transform");
                if (!label.has_value() || !sensor_id.has_value() || !transform.has_value())
                {
                    continue;
                }
                const auto sensor = sensors.find(integer(*sensor_id, "camera sensor_id"));
                if (sensor == sensors.end())
                {
                    throw std::runtime_error("Metashape camera references unknown sensor_id: " + *sensor_id);
                }
                const auto values = numericLine(*transform, 16);
                ImportedCamera camera;
                camera.imageName = *label;
                camera.calibration = sensor->second.calibration;
                camera.acquisition.rollingShutterMode = sensor->second.rollingShutterMode;
                camera.cameraToWorldRotation = {
                    values[0], values[1], values[2], values[4], values[5], values[6], values[8], values[9], values[10]};
                camera.center = {values[3], values[7], values[11]};
                cameras.push_back(std::move(camera));
            }
            if (cameras.empty())
            {
                throw std::runtime_error("Metashape document has no camera records with transforms");
            }
            return Result<std::vector<ImportedCamera>>::success(std::move(cameras));
        }
        catch (const std::exception& error)
        {
            return Result<std::vector<ImportedCamera>>::failure(
                CameraErrorCode::ParseFailure, error.what(), "Metashape XML");
        }
    }

} // namespace placamera
