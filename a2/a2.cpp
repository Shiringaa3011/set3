#include <iostream>
#include <vector>
#include <cstdint>
#include <random>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <string>

void sortH(std::vector<int64_t>& arr, int64_t left, int64_t right, size_t threshold);
void sortM(std::vector<int64_t>& arr, int64_t left, int64_t right);

class ArrayGenerator {
 private:
  std::random_device rd;
 public:
  ArrayGenerator() = default;

  std::vector<int64_t> GenerateRandomArray(size_t size, int64_t min, int64_t max) {
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int64_t> dist(min, max);

    std::vector<int64_t> res(size);
    
    for (size_t i = 0; i < size; ++i) {
      res[i] = dist(gen);
    }
    
    return res;
  }

  std::vector<int64_t> GenerateAlmostSortedArray(size_t size, int64_t min, int64_t max) {
    std::vector<int64_t> res = GenerateRandomArray(size, min, max);
    std::sort(res.begin(), res.end());
    int num_swaps = std::max(static_cast<int>(std::log2(size) * 2), 2);

    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> idx_dist(0, size - 1);
    
    for (int i = 0; i < num_swaps; ++i) {
      size_t idx1 = idx_dist(gen);
      size_t idx2 = idx_dist(gen);
      std::swap(res[idx1], res[idx2]);
    }
    
    return res;
  }

  std::vector<int64_t> GenerateReverseSortedArray(size_t size, int64_t min, int64_t max) {
    std::vector<int64_t> res = GenerateRandomArray(size, min, max);
    std::sort(res.begin(), res.end());
    std::reverse(res.begin(), res.end());
    return res;
  }
};

class SortTester {
 private:
  const size_t runs = 3;
 public:
  enum ArrayType {
    RANDOM,
    ALMOST_SORTED,
    REVERSE_SORTED
  };
    
  enum AlgorithmType {
    MERGE,
    HYBRID
  };

  double TestSortingTime(const std::vector<int64_t>& arr, AlgorithmType alg_t, size_t threshold = 0) {
    double time = 0.0;
    for (size_t i = 0; i < runs; ++i) {
      std::vector<int64_t> test_arr = arr;
      auto start = std::chrono::high_resolution_clock::now();
      if (alg_t == MERGE) {
        sortM(test_arr, 0, test_arr.size() - 1);
      } else {
        sortH(test_arr, 0, test_arr.size() - 1, threshold);
      }
      auto elapsed = std::chrono::high_resolution_clock::now()- start;
      int64_t msec = std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
      time += msec;
    }
    return time / runs;
  }
};

//insertion sort
void insertionSort(std::vector<int64_t>& arr, int64_t left, int64_t right) {
  for (int64_t i = left + 1; i <= right; ++i) {
    int64_t key = arr[i];
    int j = i - 1;

    while (j >= left && arr[j] > key) {
      arr[j + 1] = arr[j];
      j = j - 1;
    }
    arr[j + 1] = key;
  }
}

//merge sort
void sortMerge(std::vector<int64_t>& arr, int64_t left, int64_t mid, int64_t right) {
  std::vector<int64_t> sorted_part(right - left + 1);
  int64_t i = left;
  int64_t j = mid + 1;
  for (int64_t ind = 0; ind < (right - left + 1); ++ind) {
    if (i > mid || (j <= right && arr[j] < arr[i])) {
      sorted_part[ind] = arr[j];
      ++j;
    } else {
      sorted_part[ind] = arr[i];
      ++i;
    }
  }
  for (int64_t ind = left; ind <= right; ++ind) {
    arr[ind] = sorted_part[ind - left];
  }
}

//sort
void sortM(std::vector<int64_t>& arr, int64_t left, int64_t right) {
  if (left < right) {
    int64_t mid = left + (right - left) / 2;

    sortM(arr, left, mid);
    sortM(arr, mid + 1, right);

    sortMerge(arr, left, mid, right);
  }
}

void sortH(std::vector<int64_t>& arr, int64_t left, int64_t right, size_t threshold) {
  if ((right - left + 1)) {
    insertionSort(arr, left, right);
    return;
  }
  if (left < right) {
    int64_t mid = left + (right - left) / 2;

    sortH(arr, left, mid, threshold);
    sortH(arr, mid + 1, right, threshold);

    sortMerge(arr, left, mid, right);
  }
}

void SaveToCSV(const std::string& filename, const std::vector<std::vector<std::string>>& data) {
  std::ofstream file(filename);
    
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << filename << '\n';
    return;
  }
    
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
  const size_t max_size = 50000;
  const size_t step = 100;
  const size_t start_size = 40100;
  const int64_t min_val = 0;
  const int64_t max_val = 10000;
    
  ArrayGenerator generator;
  SortTester tester;

  std::vector<int64_t> rnd_arr = generator.GenerateRandomArray(max_size, min_val, max_val);
  std::vector<int64_t> almost_arr = generator.GenerateAlmostSortedArray(max_size, min_val, max_val);
  std::vector<int64_t> reverse_arr = generator.GenerateReverseSortedArray(max_size, min_val, max_val);

  std::vector<size_t> thresholds = {10, 50};

  std::vector<std::vector<std::string>> csv_data;
  std::vector<std::string> header = {"ArraySize", "ArrayType", "Algorithm", "Threshold", "TimeMs"};
  csv_data.push_back(header);

  for (size_t size = start_size; size <= max_size; size += step) {

    std::vector<int64_t> rnd_sub(rnd_arr.begin(), rnd_arr.begin() + size);
    std::vector<int64_t> almost_sub(almost_arr.begin(), almost_arr.begin() + size);
    std::vector<int64_t> reverse_sub(reverse_arr.begin(), reverse_arr.begin() + size);
        
    double time_rnd = tester.TestSortingTime(rnd_sub, SortTester::MERGE);
    double time_almost = tester.TestSortingTime(almost_sub, SortTester::MERGE);
    double time_reverse = tester.TestSortingTime(reverse_sub, SortTester::MERGE);
        
    csv_data.push_back({std::to_string(size), "Random", "StandardMerge", "0", std::to_string(time_rnd)});
    csv_data.push_back({std::to_string(size), "AlmostSorted", "StandardMerge", "0", std::to_string(time_almost)});
    csv_data.push_back({std::to_string(size), "ReverseSorted", "StandardMerge", "0", std::to_string(time_reverse)});
        
    for (size_t threshold : thresholds) {
      time_rnd = tester.TestSortingTime(rnd_sub, SortTester::HYBRID, threshold);
      time_almost = tester.TestSortingTime(almost_sub, SortTester::HYBRID, threshold);
      time_reverse = tester.TestSortingTime(reverse_sub, SortTester::HYBRID, threshold);
            
      csv_data.push_back({std::to_string(size), "Random", "Hybrid", std::to_string(threshold), std::to_string(time_rnd)});
      csv_data.push_back({std::to_string(size), "AlmostSorted", "Hybrid", std::to_string(threshold), std::to_string(time_almost)});
      csv_data.push_back({std::to_string(size), "ReverseSorted", "Hybrid", std::to_string(threshold), std::to_string(time_reverse)});
    }
  }

  SaveToCSV("sorting_results5.csv", csv_data);
  return 0;
}
