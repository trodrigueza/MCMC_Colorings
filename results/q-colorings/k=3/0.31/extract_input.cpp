#include <bits/stdc++.h>
using namespace std;

vector<pair<int, long long>> act = {{2,2},{3,246},{4,9612},{5,142820},{6,1166910},{7,6464682},{8,27350456},{9,95004072},{10,283982490},{11,754324670},{12,1821684612},{13,4067709516},{14,8506024982},{15,16822697010}};
int main() {
  string s;
  string steps;
  while (true) {
    getline(cin, s);
    if (s[0] == 'f') return 0;
    if (s[0] == 's') {
      steps = string(s.begin()+23, s.end());
    }
    if (s[0] == 'E') {
      string q = "";
      q += s[20];
      if (s[21] != '-') q += s[21];
      string est = "";
      for (int i = 49; i < s.size(); i++) {
        est += s[i];
      }

      cout << "3 " << q << " " << steps << " " << est << " " << act[stoi(q)-2].second << " " << setprecision(4) << 100*abs((long double)(act[stoi(q)-2].second-stoll(est)))/(long double)(act[stoi(q)-2].second) << " ";
      getline(cin, s);
      string time = "";
      for (int i = 12; i < s.size(); i++) {
        if (s[i] == ' ') break;
        time += s[i];
      }
      cout << time << "\n";
    }
  }
}
