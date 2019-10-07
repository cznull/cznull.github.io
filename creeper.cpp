#include <string>
#include <iostream>

int isnuorch(char x) {
	return ('0' <= x && x <= '9') || ('A' <= x && x <= 'Z') || ('a' <= x && x <= 'z');
}

void reduce(std::string &s) {
	for (int j = 0; j < s.length(); j++) {
		if (!isnuorch(s[j])) {
			s.erase(j, 1);
			j--;
		}
	}
	for (int j = 0; j < s.length(); j++) {
		if ('A' <= s[j] && s[j] <= 'Z') {
			s[j] = s[j] - 'A' + 'a';
		}
	}
}

int main()
{
	std::vector<std::string> ir = {
		"Creeper?",
		"Awww man!"
	};
	std::cout << ir[0];
  for(int j=0;j<
	for (int i = 1; i < ir.size(); i++) {
		std::string s;
		std::cin >> s;
		reduce(s);
		reduce(ir[i]);
		if (s == ir[i]) {
			i++;
			if (i < ir.size()) {
				std::cout << ir[i] << "\n";
			}
		}
		else {
			printf("please stop your kongxia behavior");
		}
	}
	return 0;
}
