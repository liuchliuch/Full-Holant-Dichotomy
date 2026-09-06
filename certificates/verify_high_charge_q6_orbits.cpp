#include <algorithm>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

constexpr int ZERO = -64;

struct SignedDsu {
  std::array<int, 20> parent{}, size{}, parity{};
  std::array<bool, 20> inconsistent{};

  SignedDsu() {
    for (int i = 0; i < 20; ++i) {
      parent[i] = i;
      size[i] = 1;
    }
  }

  std::pair<int, int> find(int v) {
    if (parent[v] == v) return {v, 0};
    const auto [root, up] = find(parent[v]);
    return {root, parity[v] ^ up};
  }

  // Add the equation q_u = -q_v.
  void add_negative_edge(int u, int v) {
    auto [ru, xu] = find(u);
    auto [rv, xv] = find(v);
    if (ru == rv) {
      if ((xu ^ xv) != 1) inconsistent[ru] = true;
      return;
    }
    if (size[ru] < size[rv]) {
      std::swap(ru, rv);
      std::swap(xu, xv);
    }
    parent[rv] = ru;
    parity[rv] = xu ^ xv ^ 1;
    size[ru] += size[rv];
    inconsistent[ru] = inconsistent[ru] || inconsistent[rv];
  }
};

struct Key {
  std::array<signed char, 20> value{};
  bool operator==(const Key& other) const { return value == other.value; }
  bool operator<(const Key& other) const { return value < other.value; }
};

struct KeyHash {
  std::size_t operator()(const Key& key) const {
    std::size_t h = 0;
    for (const signed char x : key.value) {
      h = h * 131 + static_cast<unsigned char>(x + 64);
    }
    return h;
  }
};

Key normalize(const std::array<int, 20>& raw) {
  std::map<int, std::vector<std::pair<int, int>>> components;
  for (int i = 0; i < 20; ++i) {
    if (raw[i] != 0) {
      components[std::abs(raw[i])].push_back({i, raw[i] > 0 ? 1 : -1});
    }
  }

  Key result;
  result.value.fill(ZERO);
  for (auto& [unused, members] : components) {
    (void)unused;
    std::sort(members.begin(), members.end());
    const int anchor = members.front().first;
    const int anchor_sign = members.front().second;
    for (const auto [index, sign] : members) {
      result.value[index] = static_cast<signed char>(
          (anchor + 1) * anchor_sign * sign);
    }
  }
  return result;
}

std::string encode(const Key& key) {
  std::ostringstream out;
  for (int i = 0; i < 20; ++i) {
    if (i != 0) out << ',';
    out << static_cast<int>(key.value[i]);
  }
  return out.str();
}

}  // namespace

int main(int argc, char** argv) {
  std::string orbit_output_path;
  if (argc == 3 && std::string(argv[1]) == "--emit-orbits") {
    orbit_output_path = argv[2];
  } else if (argc != 1) {
    std::cerr << "usage: verify_high_charge_q6_orbits "
                 "[--emit-orbits OUTPUT]\n";
    return 2;
  }

  std::vector<int> triples;
  for (int mask = 0; mask < 64; ++mask) {
    if (__builtin_popcount(static_cast<unsigned>(mask)) == 3) {
      triples.push_back(mask);
    }
  }
  std::map<int, int> triple_index;
  for (int i = 0; i < 20; ++i) triple_index[triples[i]] = i;

  using Edge = std::pair<int, int>;
  using Option = std::array<Edge, 2>;
  std::vector<std::array<Option, 3>> zero_pair_options;

  for (int i = 0; i < 6; ++i) {
    for (int j = i + 1; j < 6; ++j) {
      std::vector<int> external;
      for (int k = 0; k < 6; ++k) {
        if (k != i && k != j) external.push_back(k);
      }

      std::set<std::array<int, 4>> seen;
      std::vector<Option> options;
      for (int a = 0; a < 4; ++a) {
        for (int b = a + 1; b < 4; ++b) {
          const int r = (1 << external[a]) | (1 << external[b]);
          int external_mask = 0;
          for (const int k : external) external_mask |= 1 << k;
          const int rc = external_mask ^ r;

          Edge first = {triple_index[r | (1 << i)],
                        triple_index[r | (1 << j)]};
          Edge second = {triple_index[rc | (1 << i)],
                         triple_index[rc | (1 << j)]};
          if (first.first > first.second) std::swap(first.first, first.second);
          if (second.first > second.second) {
            std::swap(second.first, second.second);
          }
          if (second < first) std::swap(first, second);

          const std::array<int, 4> signature = {
              first.first, first.second, second.first, second.second};
          if (seen.insert(signature).second) options.push_back({first, second});
        }
      }
      if (options.size() != 3) {
        std::cerr << "internal error: a loop direction has " << options.size()
                  << " zero-pair options\n";
        return 1;
      }
      zero_pair_options.push_back({options[0], options[1], options[2]});
    }
  }

  long long number_of_choices = 1;
  for (int i = 0; i < 15; ++i) number_of_choices *= 3;

  std::unordered_set<Key, KeyHash> spaces;
  for (long long code = 0; code < number_of_choices; ++code) {
    long long word = code;
    SignedDsu dsu;
    for (int direction = 0; direction < 15; ++direction) {
      const Option& option = zero_pair_options[direction][word % 3];
      word /= 3;
      dsu.add_negative_edge(option[0].first, option[0].second);
      dsu.add_negative_edge(option[1].first, option[1].second);
    }

    std::map<int, std::vector<std::pair<int, int>>> components;
    for (int v = 0; v < 20; ++v) {
      const auto [root, parity] = dsu.find(v);
      components[root].push_back({v, parity});
    }

    std::array<int, 20> raw{};
    for (const auto& [root, members] : components) {
      if (dsu.inconsistent[root]) continue;
      const auto anchor = *std::min_element(members.begin(), members.end());
      for (const auto [index, parity] : members) {
        raw[index] = (anchor.first + 1) *
                     (parity == anchor.second ? 1 : -1);
      }
    }
    spaces.insert(normalize(raw));
  }

  std::vector<std::array<int, 6>> permutations;
  std::array<int, 6> permutation = {0, 1, 2, 3, 4, 5};
  do {
    permutations.push_back(permutation);
  } while (std::next_permutation(permutation.begin(), permutation.end()));

  std::set<Key> orbit_representatives;
  for (const Key& key : spaces) {
    Key best = key;
    for (const auto& p : permutations) {
      for (int complement = 0; complement < 2; ++complement) {
        std::array<int, 20> transformed{};
        for (int i = 0; i < 20; ++i) {
          if (key.value[i] == ZERO) continue;
          int mask = 0;
          for (int coordinate = 0; coordinate < 6; ++coordinate) {
            if ((triples[i] >> coordinate) & 1) mask |= 1 << p[coordinate];
          }
          if (complement) mask ^= 63;
          transformed[triple_index[mask]] = key.value[i];
        }
        best = std::min(best, normalize(transformed));
      }
    }
    orbit_representatives.insert(best);
  }

  // Canonical representatives in the explicit increasing-bitmask coordinate
  // order above.  The Python stage receives the computed set through the
  // --emit-orbits handshake and reindexes it by the emitted masks.
  const std::set<std::string> expected = {
      "-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64,-64",
      "-64,-64,-64,-64,-64,-64,-64,-64,-64,10,11,-64,-64,-64,-64,-64,-64,-64,-64,-64",
      "-64,-64,-64,-64,-64,-64,-64,-64,9,-9,11,-11,-64,-64,-64,-64,-64,-64,-64,-64",
      "-64,-64,-64,-64,-64,-64,-64,-64,9,10,11,12,-64,-64,-64,-64,-64,-64,-64,-64",
      "-64,-64,-64,-64,-64,-64,7,7,-64,-7,7,-64,-7,-7,-64,-64,-64,-64,-64,-64",
      "-64,-64,-64,-64,-64,-64,7,7,-64,-7,11,-64,-11,-11,-64,-64,-64,-64,-64,-64",
      "-64,-64,-64,-64,-64,-64,7,8,-64,-8,11,-64,-11,14,-64,-64,-64,-64,-64,-64",
      "-64,-64,-64,-64,-64,-64,7,8,-64,10,11,-64,13,14,-64,-64,-64,-64,-64,-64",
      "-64,-64,-64,-64,-64,6,-6,-6,9,-64,-64,12,13,13,-13,-64,-64,-64,-64,-64",
      "-64,-64,-64,-64,-64,6,-6,8,9,-64,-64,12,13,14,-14,-64,-64,-64,-64,-64",
      "-64,-64,-64,-64,-64,6,7,8,9,-64,-64,12,13,14,15,-64,-64,-64,-64,-64",
      "1,-1,-1,-1,-1,-1,1,-1,-1,1,-1,1,1,-1,1,1,1,1,1,-1",
      "1,-1,-1,-1,-1,-1,1,-1,-1,1,1,-1,-1,1,-1,-1,-1,-1,-1,1",
      "1,-1,-1,-1,-1,-1,1,1,-1,1,-1,1,-1,-1,1,1,1,1,1,-1",
      "1,-1,-1,-1,-1,-1,1,8,-1,1,-1,1,13,-1,1,1,1,1,1,-1",
      "1,-1,-1,-1,-1,1,1,-1,1,1,1,1,-1,1,1,-1,-1,-1,-1,1",
      "1,-1,-1,-1,-1,1,1,8,-8,1,1,12,-12,1,1,-1,-1,-1,-1,1",
      "1,-1,-1,1,-1,1,-1,1,1,-1,1,-1,-1,1,-1,1,-1,1,1,-1",
      "1,-1,-1,1,-1,1,-1,1,1,10,11,-1,-1,1,-1,1,-1,1,1,-1",
      "1,-1,-1,1,-1,1,-1,8,-8,-1,1,12,-12,1,-1,1,-1,1,1,-1",
      "1,-1,-1,1,-1,1,-1,8,9,-1,1,12,13,1,-1,1,-1,1,1,-1",
      "1,-1,-1,1,-1,6,-1,1,6,-6,6,-6,-1,1,-6,1,-1,1,1,-1",
      "1,-1,-1,1,5,-5,-1,1,-5,10,-5,10,-1,1,10,-10,-1,1,1,-1",
      "1,-1,-1,1,5,-5,-1,1,-5,10,11,5,-1,1,5,-5,-1,1,1,-1",
  };

  std::set<std::string> actual;
  for (const Key& key : orbit_representatives) actual.insert(encode(key));

  if (spaces.size() != 1868 || actual != expected) {
    std::cerr << "FAIL: got " << spaces.size() << " linear spaces and "
              << actual.size() << " symmetry orbits\n";
    for (const std::string& value : actual) {
      if (!expected.count(value)) std::cerr << "unexpected: " << value << '\n';
    }
    for (const std::string& value : expected) {
      if (!actual.count(value)) std::cerr << "missing: " << value << '\n';
    }
    return 1;
  }

  if (!orbit_output_path.empty()) {
    std::ofstream out(orbit_output_path, std::ios::out | std::ios::trunc);
    if (!out) {
      std::cerr << "FAIL: cannot open orbit output " << orbit_output_path
                << '\n';
      return 1;
    }
    out << "format=high-charge-q6-orbits-v1\n";
    out << "coordinate_order=increasing-integer-bitmask;bit-k-is-port-k\n";
    out << "triple_masks=";
    for (std::size_t i = 0; i < triples.size(); ++i) {
      if (i != 0) out << ',';
      out << triples[i];
    }
    out << '\n';
    for (const std::string& value : actual) out << "orbit=" << value << '\n';
    if (!out) {
      std::cerr << "FAIL: could not finish orbit output " << orbit_output_path
                << '\n';
      return 1;
    }
    std::cout << "HANDSHAKE: wrote " << actual.size()
              << " canonical orbits to " << orbit_output_path << '\n';
  }

  std::cout << "PASS: " << number_of_choices
            << " zero-pair choices -> 1868 signed linear spaces -> 24 orbits\n";
  return 0;
}
