#include <cmath>
#include <matplot/matplot.h>
#include <chrono>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <fstream>
#include <ctime>
#include <windows.h>
#include <string>
#include <filesystem>
#include <functional>


namespace fs = std::filesystem;
//using namespace matplot;

std::vector<int> data_a;
std::vector <double> deg(7);
int M, p, m;

std::string toString(int val)
{
	std::ostringstream oss;
	oss << val;
	return oss.str();
}

void create_data_a(int n) {
	std::vector<int> data_m(1);
	data_m[0] = M;
	srand(time(0));
	int k = 1;
	int last = M;
	while (k < m) {
		if (k < m - 1) {
			data_m.push_back(std::max(1, last - M / m) + rand() % (M / m));
			last = last - M / m;
		}
		else {
			data_m.push_back(1 + rand() % (last - 1));
			last = last - M / m;
		}
		k++;
	}
	std::vector<int> data_m_res = data_m;
	for (int i = 0; i < m; i++) {
		int ind = rand() % (m - i);
		data_a.push_back(data_m[ind]);
		std::swap(data_m[ind], data_m[data_m.size() - 1]);
		data_m.pop_back();
	}
	for (int i = m; i < n; i++) {
		data_a.push_back(data_m_res[0 + rand() % m]);
	}
}

void create_tests() {
	std::string name = "input_data\\";
	fs::create_directory(name);
	int base = 2, step = 4000;
	deg[0] = 1.0;
	for (int i = 1; i < 8; i++) {
		deg[i] = deg[i - 1] * base;
	}
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			int n = 200000;
			std::ofstream out(name + std::string("\\") + std::string(toString(i) + toString(j) + ".txt"));
			while (n >= 1 && n / deg[i] >= 2 && std::ceil(n / deg[i]) / deg[j] >= 1) {
				M = std::ceil(n / deg[i]);
				m = std::ceil(M / deg[j]);
				create_data_a(n);
				out << toString(n) + "\n";
				for (int k = 0; k < n; k++) {
					out << data_a[k] << " ";
				}
				out << "\n";
				std::cout << n << " " << i << "  " << j << "\n";
				//std::cout << "\n";
				data_a.clear();
				data_a.resize(0);
				n -= step;
			}
			out.close();
		}
	}
}

void test_generator() {
	std::string name = "input_data\\";
	fs::create_directory(std::string("input_data"));
	int base = 2;
	deg[0] = 1.0;
	for (int i = 1; i < 7; i++) {
		deg[i] = deg[i - 1] * base;
	}
	for (int i = 0; i < 7; i += 3) {
		for (int j = 0; j < 7; j += 3) {
			int n = 200000;
			std::ofstream input(name + std::string("\\") + std::string(toString(i) + toString(j) + ".txt"));
			while (n >= 1 && n / deg[i] >= 2 && std::ceil(n / deg[i]) / deg[j] >= 1) {
				M = std::ceil(n / deg[i]);
				m = std::ceil(M / deg[j]);
				create_data_a(n);
				input << toString(n) + "\n";
				for (int k = 0; k < n; k++) {
					input << data_a[k] << " ";
				}
				input << "\n";
				std::cout << n << " " << i << "  " << j << "\n";
				data_a.clear();
				data_a.resize(0);
				n -= 4000;
			}
			input.close();
		}
	}
}


class Tester {
public:
	double time_test(auto function) {
		auto begin = std::chrono::steady_clock::now();
		function();
		auto end = std::chrono::steady_clock::now();
		return double(std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count());
	}
};

void delete_bad_values(std::vector<double>& a) {
	int step = 5;
	std::vector<double> med(step);
	for (int i = 0; i < a.size(); i += step) {
		step = 5;
		double rms;
		double sum = 0;
		double sumkv = 0;
		int k = 0;
		int max_ind = std::min(int(a.size()), i + step);
		step = max_ind - i;
		for (int j = i; j < max_ind; j++) {
			sum += a[j];
			med[k] = a[j];
			k++;
		}
		std::sort(med.begin(), med.end());
		sum /= step;
		for (int j = i; j < max_ind; j++) {
			sumkv += (a[j] - sum) * (a[j] - sum);
		}
		rms = std::sqrt(sumkv / step);
		for (int j = i; j < max_ind; j++) {
			if (std::abs(med[step / 2] - a[j]) / rms >= 1.1) {
				a[j] = med[step / 2];
			}
		}
	}
}

void generate_graphs(std::string& filename, std::vector<std::vector<double>> times, std::vector<double>& n_values) {
	auto f = matplot::figure(true);
	matplot::plot(n_values, times[0], "g", n_values, times[1], "r", n_values, times[2], "b");
	matplot::ylim({ 0, 0.0025 });
	matplot::ylabel("t / n");
	matplot::xlabel("n");
	save(f, filename + ".png");
}

void read_data(int n, std::ifstream& fin, std::vector<int>& data) {
	data.resize(n);
	for (int i = 0; i < n; i++) {
		fin >> data[i];
	}
}

void s1(std::vector<int>& data) {
	std::vector<int> as(200001, -1);
	for (int i = 0; i < data.size(); i++) {
		int idx = data[i];
		as[idx] = data.size() - i - 1;
	}
	int answer = std::distance(as.begin(), std::max_element(as.begin(), as.end()));
}

void s2(std::vector<int>& data, bool reserve, int n) {
	std::unordered_set<int> unique;
	if(reserve)
		unique.reserve(200000);
	int idx_unique = data.size();
	for (int i = data.size() - 1; 0 <= i; i -= 1) {
		if (!unique.contains(data[i])) {
			idx_unique = data[i];
			unique.insert(idx_unique);
		}
	}
}

void stat() {
	std::vector<std::vector<double>> times(3);
	std::vector<double> data_n;
	int n;
	bool reserve = false;
	Tester t;
	std::string gr_name = "../../../img\\\\graphic";
	std::vector<int> test_data;
	for (int i = 0; i < 7; i += 3) {
		for (int j = 0; j < 7; j += 3) {
			std::ifstream fin(std::string("input_data\\") + toString(i) + toString(j) + ".txt");
			while (fin >> n) {
				read_data(n, fin, test_data);
				auto test_s1 = std::bind(s1, test_data);
				auto test_s2 = std::bind(s2, test_data, true, 200000);
				auto test_s3 = std::bind(s2, test_data, false, 200000);
				std::vector<double> diff_times;
				for (int k = 0; k < 3; k++) {
					diff_times.push_back(t.time_test(test_s1) / n);
				}
				std::sort(diff_times.begin(), diff_times.end());
				times[0].push_back(diff_times[1]);
				diff_times.resize(0);
				for (int k = 0; k < 3; k++) {
					diff_times.push_back(t.time_test(test_s2) / n);
				}
				std::sort(diff_times.begin(), diff_times.end());
				times[1].push_back(diff_times[1]);
				diff_times.resize(0);
				for (int k = 0; k < 3; k++) {
					diff_times.push_back(t.time_test(test_s3) / n);
				}
				std::sort(diff_times.begin(), diff_times.end());
				times[2].push_back(diff_times[1]);
				diff_times.resize(0);
				data_n.push_back(n);
				test_data.resize(0);
			}
			fin.close();
			for (int k = 0; k < times.size(); k++) {
				delete_bad_values(times[k]);
			}
			gr_name += toString(i) + "_" + toString(j);
			generate_graphs(gr_name, times, data_n);
			gr_name = "../../../img\\\\graphic";
			data_n.resize(0);
			for (int k = 0; k < times.size(); k++) {
				times[k].resize(0);
			}
		}
	}
}

int main() {
	test_generator();
	stat();
	return 0;
}
