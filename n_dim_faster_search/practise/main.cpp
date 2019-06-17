#include <iostream>
#include <vector>

using namespace std;

int main() {

	vector<int> a = {2,5,3,8,2,5,7,6};

	auto iter = std::find(a.begin(), a.end(), 2);

	for (iter; iter != a.end(); iter++)
	{
		std::cout << *iter << std::endl;
	}

	system("pause");
	return 0;
}