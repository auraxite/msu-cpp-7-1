#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {

using u64 = uint64_t;
using u128 = unsigned __int128;

constexpr u64 kTop = 1ULL << 56;
constexpr int kProbBits = 12;
constexpr int kAdaptShift = 5;

class Encoder {
 public:
  void Encode(u64 cum, u64 freq, u64 total) {
    u64 r = range_ / total;
    low_ += static_cast<u128>(r * cum);
    range_ = r * freq;
    while (range_ < kTop) {
      ShiftLow();
      range_ <<= 8;
    }
  }
  void EncodeRaw(u64 value, int bits) {
    for (; bits > 16; bits -= 16) {
      Encode((value >> (bits - 16)) & 0xFFFF, 1, 1 << 16);
    }
    Encode(value & ((1ULL << bits) - 1), 1, 1ULL << bits);
  }
  void EncodeBit(uint16_t* p, int bit) {
    if (bit) {
      Encode(*p, (1 << kProbBits) - *p, 1 << kProbBits);
      *p -= *p >> kAdaptShift;
    } else {
      Encode(0, *p, 1 << kProbBits);
      *p += ((1 << kProbBits) - *p) >> kAdaptShift;
    }
  }
  std::string Finish() {
    for (int i = 0; i < 9; ++i) ShiftLow();
    return out_;
  }

 private:
  void ShiftLow() {
    if (static_cast<u64>(low_) < 0xFF00000000000000ULL || (low_ >> 64) != 0) {
      uint8_t carry = static_cast<uint8_t>(low_ >> 64);
      uint8_t temp = cache_;
      do {
        out_.push_back(static_cast<char>(temp + carry));
        temp = 0xFF;
      } while (--cache_size_ != 0);
      cache_ = static_cast<uint8_t>(static_cast<u64>(low_) >> 56);
    }
    ++cache_size_;
    low_ = static_cast<u64>(static_cast<u64>(low_) << 8);
  }

  u128 low_ = 0;
  u64 range_ = ~0ULL;
  uint8_t cache_ = 0;
  u64 cache_size_ = 1;
  std::string out_;
};

class Decoder {
 public:
  explicit Decoder(const std::string& in) : in_(in) {
    for (int i = 0; i < 9; ++i) code_ = (code_ << 8) | NextByte();
  }
  u64 Target(u64 total) {
    r_ = range_ / total;
    return std::min(code_ / r_, total - 1);
  }
  void Update(u64 cum, u64 freq) {
    code_ -= r_ * cum;
    range_ = r_ * freq;
    while (range_ < kTop) {
      code_ = (code_ << 8) | NextByte();
      range_ <<= 8;
    }
  }
  u64 DecodeRaw(int bits) {
    u64 value = 0;
    for (; bits > 0; bits -= 16) {
      int n = std::min(bits, 16);
      u64 t = Target(1ULL << n);
      Update(t, 1);
      value = (value << n) | t;
    }
    return value;
  }
  int DecodeBit(uint16_t* p) {
    if (Target(1 << kProbBits) >= *p) {
      Update(*p, (1 << kProbBits) - *p);
      *p -= *p >> kAdaptShift;
      return 1;
    }
    Update(0, *p);
    *p += ((1 << kProbBits) - *p) >> kAdaptShift;
    return 0;
  }

 private:
  u64 NextByte() {
    return pos_ < in_.size() ? static_cast<uint8_t>(in_[pos_++]) : 0;
  }

  const std::string& in_;
  size_t pos_ = 0;
  u64 range_ = ~0ULL;
  u64 code_ = 0;
  u64 r_ = 1;
};

class Fenwick {
 public:
  explicit Fenwick(const std::vector<uint32_t>& w) : tree_(w.size() + 1) {
    for (size_t i = 0; i < w.size(); ++i) tree_[i + 1] = w[i];
    for (size_t i = 1; i < tree_.size(); ++i) {
      size_t j = i + (i & -i);
      if (j < tree_.size()) tree_[j] += tree_[i];
    }
    while (high_bit_ * 2 < tree_.size()) high_bit_ *= 2;
  }
  void Add(size_t i, int64_t delta) {
    for (++i; i < tree_.size(); i += i & -i) tree_[i] += delta;
  }
  u64 Prefix(size_t n) const {
    u64 s = 0;
    for (; n > 0; n -= n & -n) s += tree_[n];
    return s;
  }
  size_t Find(u64 x, u64* rest) const {
    size_t pos = 0;
    for (size_t step = high_bit_; step > 0; step >>= 1) {
      if (pos + step < tree_.size() && tree_[pos + step] <= x) {
        pos += step;
        x -= tree_[pos];
      }
    }
    *rest = x;
    return pos;
  }

 private:
  std::vector<u64> tree_;
  size_t high_bit_ = 1;
};

struct Models {
  uint16_t gap[64];
  uint16_t len[33];
  uint16_t mant[33][33];
  uint16_t loop;
  Models() {
    std::fill_n(&gap[0], 64, 1 << (kProbBits - 1));
    std::fill_n(&len[0], 33, 1 << (kProbBits - 1));
    std::fill_n(&mant[0][0], 33 * 33, 1 << (kProbBits - 1));
    loop = 1 << (kProbBits - 1);
  }
};

int RiceShift(u64 n) {
  int k = 0;
  while (k < 31 && (n << (k + 1)) <= (1ULL << 32) * 7 / 10) ++k;
  return k;
}

void EncodeGamma(Encoder* enc, Models* m, u64 n) {
  int bits = 64 - __builtin_clzll(n);
  for (int j = 1; j < 33; ++j) {
    enc->EncodeBit(&m->len[j], bits > j);
    if (bits == j) break;
  }
  for (int j = bits - 2; j >= 0; --j) {
    enc->EncodeBit(&m->mant[bits][j], (n >> j) & 1);
  }
}

u64 DecodeGamma(Decoder* dec, Models* m) {
  int bits = 1;
  while (bits < 33 && dec->DecodeBit(&m->len[bits])) ++bits;
  u64 n = 1;
  for (int j = bits - 2; j >= 0; --j) {
    n = (n << 1) | dec->DecodeBit(&m->mant[bits][j]);
  }
  return n;
}

bool ReadFile(const char* path, std::string* data) {
  FILE* f = std::fopen(path, "rb");
  if (!f) return false;
  std::fseek(f, 0, SEEK_END);
  data->resize(std::ftell(f));
  std::fseek(f, 0, SEEK_SET);
  size_t got = std::fread(&(*data)[0], 1, data->size(), f);
  std::fclose(f);
  return got == data->size();
}

bool WriteFile(const char* path, const std::string& data) {
  FILE* f = std::fopen(path, "wb");
  if (!f) return false;
  size_t put = std::fwrite(data.data(), 1, data.size(), f);
  std::fclose(f);
  return put == data.size();
}

std::string Serialize(const std::string& text) {
  std::vector<uint32_t> raw;
  u64 cur = 0;
  bool in_num = false;
  for (char c : text) {
    if (c >= '0' && c <= '9') {
      cur = cur * 10 + (c - '0');
      in_num = true;
    } else if (in_num) {
      raw.push_back(static_cast<uint32_t>(cur));
      cur = 0;
      in_num = false;
    }
  }
  if (in_num) raw.push_back(static_cast<uint32_t>(cur));
  size_t edges = raw.size() / 3;

  std::vector<uint32_t> ids;
  for (size_t i = 0; i < edges; ++i) {
    ids.push_back(raw[3 * i]);
    ids.push_back(raw[3 * i + 1]);
  }
  std::sort(ids.begin(), ids.end());
  ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
  size_t n = ids.size();
  auto index = [&](uint32_t id) {
    auto it = std::lower_bound(ids.begin(), ids.end(), id);
    return static_cast<uint32_t>(it - ids.begin());
  };

  std::vector<uint32_t> deg(n, 0), start(n + 1, 0);
  std::vector<int> loop_weight(n, -1);
  std::vector<u64> packed;  // (a << 40) | (b << 8) | w with a < b
  for (size_t i = 0; i < edges; ++i) {
    uint32_t a = index(raw[3 * i]);
    uint32_t b = index(raw[3 * i + 1]);
    uint32_t w = raw[3 * i + 2];
    if (a == b) {
      loop_weight[a] = static_cast<int>(w);
      continue;
    }
    if (a > b) std::swap(a, b);
    ++deg[a];
    ++deg[b];
    packed.push_back((static_cast<u64>(a) << 40) | (static_cast<u64>(b) << 8) |
                     w);
  }
  std::vector<uint32_t>().swap(raw);
  std::sort(packed.begin(), packed.end());

  Encoder enc;
  Models m;
  enc.EncodeRaw(n, 32);
  if (n == 0) return enc.Finish();
  int k = RiceShift(n);
  enc.EncodeRaw(ids[0], 32);
  for (size_t i = 1; i < n; ++i) {
    u64 g = ids[i] - ids[i - 1] - 1;
    u64 q = g >> k;
    for (u64 j = 0; j <= q; ++j) {
      enc.EncodeBit(&m.gap[std::min<u64>(j, 63)], j < q);
    }
    if (k > 0) enc.EncodeRaw(g & ((1ULL << k) - 1), k);
  }
  for (size_t u = 0; u < n; ++u) {
    EncodeGamma(&enc, &m, deg[u] + 1);
    enc.EncodeBit(&m.loop, loop_weight[u] >= 0);
    if (loop_weight[u] >= 0) enc.EncodeRaw(loop_weight[u], 8);
  }

  Fenwick fw(deg);
  u64 sum = packed.size() * 2;
  size_t e = 0;
  for (size_t u = 0; u < n; ++u) {
    fw.Add(u, -static_cast<int64_t>(deg[u]));
    sum -= deg[u];
    size_t prev = u;
    for (; e < packed.size() && (packed[e] >> 40) == u; ++e) {
      size_t v = (packed[e] >> 8) & 0xFFFFFFFF;
      u64 base = fw.Prefix(prev + 1);
      enc.Encode(fw.Prefix(v) - base, deg[v], sum - base);
      enc.EncodeRaw(packed[e] & 0xFF, 8);
      fw.Add(v, -1);
      --deg[v];
      --sum;
      prev = v;
    }
  }
  return enc.Finish();
}

void AppendLine(std::string* out, uint32_t a, uint32_t b, uint32_t w) {
  char buf[40];
  int len = std::snprintf(buf, sizeof(buf), "%u\t%u\t%u\n", a, b, w);
  out->append(buf, len);
}

std::string Deserialize(const std::string& bin) {
  Decoder dec(bin);
  Models m;
  std::string out;
  size_t n = dec.DecodeRaw(32);
  if (n == 0) return out;
  int k = RiceShift(n);
  std::vector<uint32_t> ids(n), deg(n);
  ids[0] = dec.DecodeRaw(32);
  for (size_t i = 1; i < n; ++i) {
    u64 q = 0;
    while (dec.DecodeBit(&m.gap[std::min<u64>(q, 63)])) ++q;
    u64 g = (q << k) | (k > 0 ? dec.DecodeRaw(k) : 0);
    ids[i] = static_cast<uint32_t>(ids[i - 1] + g + 1);
  }
  for (size_t u = 0; u < n; ++u) {
    deg[u] = DecodeGamma(&dec, &m) - 1;
    if (dec.DecodeBit(&m.loop)) {
      AppendLine(&out, ids[u], ids[u], dec.DecodeRaw(8));
    }
  }

  Fenwick fw(deg);
  u64 sum = 0;
  for (uint32_t d : deg) sum += d;
  for (size_t u = 0; u < n; ++u) {
    uint32_t r = deg[u];
    fw.Add(u, -static_cast<int64_t>(r));
    sum -= r;
    size_t prev = u;
    for (uint32_t i = 0; i < r; ++i) {
      u64 base = fw.Prefix(prev + 1);
      u64 rest;
      u64 target = dec.Target(sum - base);
      size_t v = fw.Find(target + base, &rest);
      dec.Update(target - rest, deg[v]);
      AppendLine(&out, ids[u], ids[v], dec.DecodeRaw(8));
      fw.Add(v, -1);
      --deg[v];
      --sum;
      prev = v;
    }
  }
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  const char* mode = nullptr;
  const char* input = nullptr;
  const char* output = nullptr;
  for (int i = 1; i < argc; ++i) {
    if (!std::strcmp(argv[i], "-s") || !std::strcmp(argv[i], "-d")) {
      mode = argv[i];
    } else if (!std::strcmp(argv[i], "-i") && i + 1 < argc) {
      input = argv[++i];
    } else if (!std::strcmp(argv[i], "-o") && i + 1 < argc) {
      output = argv[++i];
    }
  }
  if (!mode || !input || !output) {
    std::fprintf(stderr, "usage: %s -s|-d -i input -o output\n", argv[0]);
    return 1;
  }
  std::string data;
  if (!ReadFile(input, &data)) {
    std::fprintf(stderr, "cannot read %s\n", input);
    return 1;
  }
  std::string result = mode[1] == 's' ? Serialize(data) : Deserialize(data);
  if (!WriteFile(output, result)) {
    std::fprintf(stderr, "cannot write %s\n", output);
    return 1;
  }
  return 0;
}
