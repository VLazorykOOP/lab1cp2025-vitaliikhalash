#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>

std::vector<std::pair<double, double>> Ux_data;
std::vector<std::pair<double, double>> Tx_data;
std::vector<std::pair<std::string, double>> Ctext_data;

class FileOpenException
{
    std::string filename;

public:
    explicit FileOpenException(const std::string &fname) : filename(fname) {}

    void Message() const
    {
        std::cout << "File could not be opened: " << filename << std::endl;
    }
};

std::vector<std::pair<double, double>> loadData(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw FileOpenException(filename);
    }

    std::vector<std::pair<double, double>> data;
    std::string line;

    std::getline(file, line);

    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string xStr;
        std::string uStr;
        if (std::getline(ss, xStr, ',') && std::getline(ss, uStr))
        {
            double x = std::stod(xStr);
            double u = std::stod(uStr);
            data.emplace_back(x, u);
        }
    }

    return data;
}

std::vector<std::pair<std::string, double>> loadData2(const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw FileOpenException(filename);
    }

    std::vector<std::pair<std::string, double>> data;
    std::string line;

    std::getline(file, line);

    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::string xStr;
        std::string uStr;
        if (std::getline(ss, xStr, ',') && std::getline(ss, uStr))
        {
            double u = std::stod(uStr);
            data.emplace_back(xStr, u);
        }
    }

    return data;
}

double get_Ux(double x)
{
    for (size_t i = 0; i < Ux_data.size(); ++i)
    {
        if (Ux_data[i].first == x)
            return Ux_data[i].second;
        if (Ux_data[i].first < x && x < Ux_data[i + 1].first)
            return Ux_data[i].second + (Ux_data[i + 1].second - Ux_data[i + 1].second) * (x - Ux_data[i].first) / (Ux_data[i + 1].first - x);
    }

    throw std::out_of_range("X is outside the data range");
}

double get_Tx(double x)
{
    for (size_t i = 0; i < Tx_data.size(); ++i)
    {
        if (Tx_data[i].first == x)
            return Tx_data[i].second;
        if (Ux_data[i].first < x && x < Ux_data[i + 1].first)
            return Tx_data[i].second + (Tx_data[i + 1].second - Tx_data[i + 1].second) * (x - Tx_data[i].first) / (Tx_data[i + 1].first - x);
    }

    throw std::out_of_range("X is outside the data range");
}

double get_Ctext(const std::string_view &text)
{
    for (const auto &[key, value] : Ctext_data)
    {
        if (key == text)
        {
            return value;
        }
    }

    return 0;
}

double r1_U(double x) { return get_Ux(x); }
double r1_T(double x) { return get_Tx(x); }
double r1_Qqn(double x, double y, double z) { return x / r1_U(x) + y * r1_T(y) - r1_U(z) * r1_T(z); }
double r1_Qnk(double x, double y) { return r1_Qqn(x, y, x + y) - r1_Qqn(y, x, x - y); }
double r1_Rnk(double x, double y) { return x * r1_Qnk(x, y) + y * r1_Qnk(y, x); }
double r1_func(double x, double y, double z) { return r1_Rnk(x, y) + r1_Rnk(y, z) * r1_Rnk(x, y); }

double r2_U(double x) { return atan(asin(cos(x * 3))); }
double r2_T(double x) { return atan(acos(sin(x * 2))); }
double r2_Qqn(double x, double y, double z) { return x / r2_U(x) + y * r2_T(y) - r2_U(z) * r2_T(z); }
double r2_Qnk(double x, double y) { return 1.1 * r2_Qqn(x, y, x + y) - 0.9 * r2_Qqn(y, x, x - y); }
double r2_Rnk(double x, double y) { return x * r2_Qnk(x, y) + y * r2_Qnk(y, x) - 0.03 * r2_Qnk(x, y) * r2_Qnk(y, x); }
double r2_func(double x, double y, double z) { return r2_Rnk(x, y) + r2_Rnk(y, z) * r2_Rnk(x, y); }

double r3_U(double x) { return get_Ux(x); }
double r3_T(double x) { return get_Tx(x); }
double r3_Qqn(double x, double y, double z) { return x / r2_U(x) + y * r2_T(y) - 0.9 * r2_U(z) * r2_T(z); }
double r3_Qnk(double x, double y) { return 1.3 * r2_Qqn(x, y, x) - r2_Qqn(y, x, x); }
double r3_func(double x, double y) { return 1.75 * x * r3_Qnk(x, y) + 1.25 * y * r3_Qnk(y, x) - 1.5 * r3_Qnk(x, y) * r3_Qnk(y, x); }

double k_Ctext(double x, const std::string &text)
{
    if (text.empty())
        return get_Ctext("set") + get_Ctext("get") - x;
    if (x > 0)
        return get_Ctext(text) + x;
    if (x <= 0)
        return get_Ctext("set") + get_Ctext(text);
    return 0;
}

double k_Max(double x, double y, double z, double u)
{
    double max = x;

    if (y > max)
        max = y;
    if (z > max)
        max = z;
    if (u > max)
        max = u;

    return max;
}

double k_Rtext(double x, double y, double z, const std::string &text)
{
    return k_Ctext(k_Max(x, y, x + z, y + z), text);
}

double r_func(int alg, double x, double y, double z)
{
    switch (alg)
    {
    case 1:
        return r1_func(x, y, z);
    case 2:
        return r2_func(x, y, z);
    case 3:
        return r3_func(x, y);
    default:
        throw std::invalid_argument("Invalid algorithm");
    }
}

double Variant(double r, double k)
{
    return 0.8973 * r + 0.1027 * k;
}

int main()
{
    double x;
    double y;
    double z;
    double r;
    double k;
    std::string text;
    std::string dat1_path;
    std::string dat2_path;
    std::string dat3_path;

    std::cout << "Input path to dat1.dat: ";
    std::cin >> dat1_path;
    std::cout << "Input path to dat2.dat: ";
    std::cin >> dat2_path;
    std::cout << "Input path to dat3.dat: ";
    std::cin >> dat3_path;
    std::cout << "Input text: ";
    std::cin >> text;
    std::cout << "Input x: ";
    std::cin >> x;
    std::cout << "Input y: ";
    std::cin >> y;
    std::cout << "Input z: ";
    std::cin >> z;

    try
    {
        Ux_data = loadData(dat1_path);

        if (fabs(x) <= 5)
        {
            r = r_func(2, x, y, z);
        }
    }
    catch (const FileOpenException &e)
    {
        e.Message();
        r = r_func(2, x, y, z);
    }
    try
    {
        Tx_data = loadData(dat2_path);

        if (fabs(x) <= 10)
        {
            r = r_func(2, x, y, z);
        }
    }
    catch (const FileOpenException &e)
    {
        e.Message();
        r = r_func(3, x, y, z);
    }
    try
    {
        Ctext_data = loadData2(dat3_path);
        k = k_Rtext(x, y, z, text);
    }
    catch (const FileOpenException &e)
    {
        e.Message();
        std::cerr << "Error while opening file. Exiting program" << std::endl;
        return 1;
    }

    double variant = Variant(r, k);

    std::cout << "r = " << r << std::endl;
    std::cout << "k = " << k << std::endl;
    std::cout << "Variant(r, k) = " << variant << std::endl;
}
