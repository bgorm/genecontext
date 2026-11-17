#pragma once
#include <htslib/hts.h>
#include <htslib/kstring.h>
#include <algorithm>
#include <charconv>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace gc {
namespace fs = std::filesystem;
inline void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
inline std::vector<std::string> split(std::string_view text, char separator) {
    std::vector<std::string> result;
    size_t start = 0;
    do {
        size_t end = text.find(separator, start);
        result.emplace_back(text.substr(start, end == text.npos ? end : end-start));
        if (end == text.npos) break;
        start = end+1;
    } while (true);
    return result;
}
inline std::string trim(std::string_view value) {
    const auto a = value.find_first_not_of(" \t\r\n");
    if (a == value.npos) return {};
    return std::string(value.substr(a, value.find_last_not_of(" \t\r\n")-a+1));
}
inline int64_t integer(std::string_view value, int64_t minimum, int64_t maximum,
                       const std::string& description) {
    int64_t result = 0;
    auto parsed = std::from_chars(value.data(), value.data()+value.size(), result);
    require(parsed.ec == std::errc{} && parsed.ptr == value.data()+value.size() &&
            result >= minimum && result <= maximum, "Invalid " + description + ": " + std::string(value));
    return result;
}
inline int hex_value(char c) {
    if (c >= '0' && c <= '9') return c-'0';
    if (c >= 'a' && c <= 'f') return c-'a'+10;
    if (c >= 'A' && c <= 'F') return c-'A'+10;
    return -1;
}
inline std::string unescape(std::string_view value) {
    std::string result;
    for (size_t i = 0; i < value.size(); ++i) {
        if (value[i] != '%') { result += value[i]; continue; }
        require(i+2 < value.size() && hex_value(value[i+1]) >= 0 && hex_value(value[i+2]) >= 0,
                "Malformed percent escape in annotation attribute");
        char c = static_cast<char>(16*hex_value(value[i+1])+hex_value(value[i+2]));
        require(c != '\0', "NUL in annotation attribute");
        result += c; i += 2;
    }
    return result;
}
// Escape structural delimiters in INFO lists, maps and the internal TSV.
inline std::string escape(std::string_view value, bool map = false) {
    static constexpr char hex[] = "0123456789ABCDEF";
    std::string result;
    for (unsigned char c : value) {
        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') ||
            (c >= '0' && c <= '9') || c == '_' || c == '-' || c == ':' ||
            (c == '.' && !map && value != ".")) result += static_cast<char>(c);
        else { result += '%'; result += hex[c >> 4]; result += hex[c & 15]; }
    }
    return result;
}
template <typename Range>
inline std::string join(const Range& values, const std::string& delimiter = ",") {
    std::string result;
    for (const auto& value : values) {
        if (!result.empty()) result += delimiter;
        result += value;
    }
    return result;
}

template <typename T>
inline void sort_unique(std::vector<T>& values) {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
}

struct KString {
    kstring_t value{0,0,nullptr};
    KString() = default;
    KString(const KString&) = delete;
    KString& operator=(const KString&) = delete;
    ~KString() { std::free(value.s); }
    std::string_view view() const { return {value.s ? value.s : "", value.l}; }
};
struct HtsCloser { void operator()(htsFile* p) const { if (p) hts_close(p); } };
using HtsFile = std::unique_ptr<htsFile,HtsCloser>;
}
