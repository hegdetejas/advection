#include <iostream>
#include <iterator>
#include <string>
#include <fstream>

template<typename InputIterator1, typename InputIterator2>
bool
range_equal(InputIterator1 first1, InputIterator1 last1,
        InputIterator2 first2, InputIterator2 last2)
{
    while(first1 != last1 && first2 != last2)
    {
        if(*first1 != *first2) return false;
        ++first1;
        ++first2;
    }
    return (first1 == last1) && (first2 == last2);
}

bool compare_files(const std::string& filename1, const std::string& filename2)
{
    std::ifstream file1(filename1);
    std::ifstream file2(filename2);

    std::istreambuf_iterator<char> begin1(file1);
    std::istreambuf_iterator<char> begin2(file2);

    std::istreambuf_iterator<char> end;

    return range_equal(begin1, end, begin2, end);
}

int main() {
  bool same_0 = compare_files("out_0_base.txt", "out_0_boost.txt");
  bool same_NTby2 = compare_files("out_NTby2_base.txt", "out_NTby2_boost.txt");
  bool same_NT = compare_files("out_NT_base.txt", "out_NT_boost.txt");
  if (same_0 && same_NTby2 && same_NT) {
    std::cout << "files are bitwise the same" << std::endl;
  } else {
    std::cout << "files are not the same" << std::endl;
  }

  return 0;
}
