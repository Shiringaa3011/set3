#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <fstream>

class Circle {
 public:
  double x = 0;
  double y = 0;
  double r = 0;

  Circle(double new_x, double new_y, double new_r) : x(new_x), y(new_y), r(new_r) {}
  Circle() {}
};

bool pointInAllCircles(double x, double y, const Circle& c1, const Circle& c2, const Circle& c3) {
  bool in1 = ((x - c1.x) * (x - c1.x) + (y - c1.y) * (y - c1.y) - c1.r * c1.r) <= 1e-9;
  bool in2 = ((x - c2.x) * (x - c2.x) + (y - c2.y) * (y - c2.y) - c2.r * c2.r) <= 1e-9;
  bool in3 = ((x - c3.x) * (x - c3.x) + (y - c3.y) * (y - c3.y) - c3.r * c3.r) <= 1e-9;
  return in1 && in2 && in3;
}

void findBoundingRectangle(const Circle& c1, const Circle& c2, const Circle& c3, double& x_min, double& x_max, double& y_min, double& y_max) {
  x_min = std::min({c1.x - c1.r, c2.x - c2.r, c3.x - c3.r});
  x_max = std::max({c1.x + c1.r, c2.x + c2.r, c3.x + c3.r});
  y_min = std::min({c1.y - c1.r, c2.y - c2.r, c3.y - c3.r});
  y_max = std::max({c1.y + c1.r, c2.y + c2.r, c3.y + c3.r});
}

double calculateArea(const Circle& c1, const Circle& c2, const Circle& c3, size_t points_all, double x_min, double x_max, double y_min, double y_max) {
  std::mt19937 gen(67);
  std::uniform_real_distribution<double> dist_x(x_min, x_max);
  std::uniform_real_distribution<double> dist_y(y_min, y_max);
  int points_in = 0;
  double s_rec = (x_max - x_min) * (y_max - y_min);
  for (int i = 0; i < points_all; ++i) {
    double x = dist_x(gen);
    double y = dist_y(gen);
    if (pointInAllCircles(x, y, c1, c2, c3)) {
      points_in++;
    }
  }
  double area = (static_cast<double>(points_in) / points_all) * s_rec;
  return area;
}

void writeToCSV(const std::string& filename, const std::vector<std::vector<double>>& data, const std::vector<std::string>& headers) {
  std::ofstream file(filename);
  
  for (size_t i = 0; i < headers.size(); ++i) {
    file << headers[i];
    if (i != headers.size() - 1) {
      file << ",";
    }
  }
  file << "\n";

  for (const auto& row : data) {
    for (size_t i = 0; i < row.size(); ++i) {
      file << row[i];
      if (i != row.size() - 1) {
        file << ",";
      }
    }
    file << "\n";
  }
  
  file.close();
}

int main() {
  Circle c1(1.0, 1.0, 1.0);
  Circle c2(1.5, 2.0, sqrt(5.0)/2.0);
  Circle c3(2.0, 1.5, sqrt(5.0)/2.0);

  double exact_area = 0.25 * M_PI + 1.25 * asin(0.8) - 1.0;
  std::cout << "Exact area: " << exact_area << std::endl;

  double wide_x_min = 0;
  double wide_x_max = 2.0 + sqrt(5.0)/2.0;
  double wide_y_min = 0;
  double wide_y_max = 2.0 + sqrt(5.0)/2.0;

  double narrow_x_min = 2.0 - sqrt(5.0)/2.0;
  double narrow_x_max = 2;
  double narrow_y_min = 2.0 - sqrt(5.0)/2.0;
  double narrow_y_max = 2;

  std::vector<std::vector<double>> experiment_data;
  
  for (int n = 100; n <= 100000; n += 500) {
    double area_wide = calculateArea(c1, c2, c3, n, wide_x_min, wide_x_max, wide_y_min, wide_y_max);
    double error_wide = std::abs(area_wide - exact_area) / exact_area * 100.0;
    double rectangle_area_wide = (wide_x_max - wide_x_min) * (wide_y_max - wide_y_min);

    double area_narrow = calculateArea(c1, c2, c3, n, narrow_x_min, narrow_x_max, narrow_y_min, narrow_y_max);
    double error_narrow = std::abs(area_narrow - exact_area) / exact_area * 100.0;
    double rectangle_area_narrow = (narrow_x_max - narrow_x_min) * (narrow_y_max - narrow_y_min);
    
    std::vector<double> row = {
      static_cast<double>(n), // Количество точек
      area_wide, // Приближенная площадь (широкая область)
      error_wide, // Относительная ошибка (%) (широкая область)
      rectangle_area_wide, // Площадь прямоугольной области (широкая)
      area_narrow, // Приближенная площадь (узкая область)
      error_narrow, // Относительная ошибка (%) (узкая область)
      rectangle_area_narrow, // Площадь прямоугольной области (узкая)
      exact_area // Точная площадь
    };
    experiment_data.push_back(row);
  }
  
  writeToCSV("monte_carlo_experiment.csv", experiment_data, 
             {"N", "AreaWide", "ErrorWidePercent", "RectangleAreaWide", 
              "AreaNarrow", "ErrorNarrowPercent", "RectangleAreaNarrow", "ExactArea"});
}
