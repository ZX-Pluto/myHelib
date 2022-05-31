#include<iostream>
#include<bits/stdc++.h>
#include<vector>
// #include "my_hexl.h"
using namespace std;
class Ntt {
public:
	uint64_t kesai;
	uint64_t q;
	uint64_t degree;
	uint64_t inv_kesai;
	uint64_t m_pow;
	vector<uint64_t> w_rom;
	Ntt(uint64_t q, uint64_t degree, uint64_t kesai):q(q), degree(degree), kesai(kesai) {
		m_pow = Log2(degree);
		inv_kesai = PowMod(kesai, (degree<<1) -1, q);
		w_rom.resize(degree);
		for (long i = 0; i < degree; ++i)
		{
			w_rom[i] = PowMod(kesai, bitreverse(i,m_pow) % q, q);
		}
	}

	Ntt(uint64_t q, uint64_t degree):q(q),degree(degree) {
		kesai = GeneratePrimitiveRoot(degree << 1, q);
		new(this) Ntt(q, degree, kesai);
	}

	Ntt(){}

	// void init(uint64_t q, uint64_t degree, uint64_t kesai){
	// 	q = q; degree = degree; kesai = kesai;
	// 	m_pow = Log2(degree);
	// 	inv_kesai = PowMod(kesai, (degree<<1) -1, q);
	// 	w_rom.resize(degree);
	// 	for (long i = 0; i < degree; ++i)
	// 	{
	// 		w_rom[i] = PowMod(kesai, bitreverse(i,m_pow) % q, q);
	// 	}
	// }

	// void init(uint64_t q, uint64_t degree){
	// 	kesai = GeneratePrimitiveRoot(degree << 1, q);
	// 	init(q, degree, kesai);
	// }

	uint64_t get_kesai() {return kesai;}
	uint64_t get_q() {return q;}
	uint64_t get_degree() {return degree;}
	uint64_t get_inv_kesai() {return inv_kesai;}
	uint64_t get_m_pow() { return m_pow;}

	uint64_t bitreverse(uint64_t n, uint64_t l)
	{
		uint64_t res = 0;
		for (long i = 0; i < l; ++i)
		{
			res = (res << 1) | (n & 1);
			n >>= 1;
		}
		return res;
	}
	vector<uint64_t> ntt(vector<uint64_t> a) {
		uint64_t n = a.size();
		uint64_t log_n = Log2(n);
		uint64_t r = 1;
		// ofstream fout("nwc_NTT_C.txt"); //文件输出流对象
		// streambuf* pOld =cout.rdbuf(fout.rdbuf());
		// int my_i = 0;
		for (long p = log_n - 1; p >= 0; --p)
		{
			int J = 1 << p;
			for (long k = 0; k < n/(J << 1); ++k)
			{
				uint64_t w = w_rom[r++];
				for (long j = 0; j < J; ++j)
				{
					uint64_t u = a[(k<<1)*J + j] % q;
					uint64_t t = MultiplyMod(a[(k<<1)*J + j + J], w, q);
					a[(k<<1)*J + j] = (u + t) % q;
					a[(k<<1)*J + j + J] = (u + q - t) % q;
				}	
			}
		}
		return a;
	}

	uint64_t op21(uint64_t a)
	{
		return (a & 1) == 0 ? (a >> 1) % q : ((a >> 1) + ((q + 1) >> 1)) % q;
	}

	vector<uint64_t> rntt(vector<uint64_t> a) {
		uint64_t n = a.size();
		uint64_t log_n = Log2(n);
		uint64_t r = w_rom.size() - 1;
		for (long i = 0; i < log_n; ++i)
		{
			uint64_t J = 1 << i;
			for (long k = 0; k < n/(J << 1); ++k)
			{
				uint64_t w = w_rom[r--];
				for (long j = 0; j < J; ++j)
				{
					uint64_t u = a[(k<<1)*J + j] % q;
					uint64_t t = a[(k<<1)*J + j + J] % q;
					a[(k<<1)*J + j] = (op21(u + t)) % q;
					a[(k<<1)*J + j + J] = MultiplyMod(op21(t + q - u), w, q);
				}	
			}
		}
		return a;
	}
	

	// typedef long long uint64_t;
	typedef __int128 int128_t;
	typedef unsigned __int128 uint128_t;

	// Returns base^exp mod modulus
	uint64_t PowMod(uint64_t base, uint64_t exp, uint64_t modulus) {
	base %= modulus;
	uint64_t result = 1;
	while (exp > 0) {
		if (exp & 1) {
		result = MultiplyMod(result, base, modulus);
		}
		base = MultiplyMod(base, base, modulus);
		exp >>= 1;
	}
	return result;
	}

	// Returns true whether root is a degree-th root of unity
	// degree must be a power of two.
	bool IsPrimitiveRoot(uint64_t root, uint64_t degree, uint64_t modulus) {
	if (root == 0) {
		return false;
	}
	if(!IsPowerOfTwo(degree)){ std::cout << degree << " not a power of 2" << std::endl; exit(0);}
	// Check if root^(degree/2) == -1 mod modulus
	return PowMod(root, degree / 2, modulus) == (modulus - 1);
	}

	/// @brief Returns whether or not num is a power of two
	bool IsPowerOfTwo(uint64_t num) { return num && !(num & (num - 1)); }

	// /// @brief Returns floor(log2(x))
	// inline uint64_t Log2(uint64_t x) { return MSB(x); }

	// inline bool IsPowerOfFour(uint64_t num) {
	//   return IsPowerOfTwo(num) && (Log2(num) % 2 == 0);
	// }

	// /// @brief Returns the maximum value that can be represented using \p bits bits
	// inline uint64_t MaximumValue(uint64_t bits) {
	//   HEXL_CHECK(bits <= 64, "MaximumValue requires bits <= 64; got " << bits);
	//   if (bits == 64) {
	//     return (std::numeric_limits<uint64_t>::max)();
	//   }
	//   return (1ULL << bits) - 1;
	// }

	// uint64_t BarrettReduce128_mod(uint64_t x, uint64_t modulus)
	// {
	
	// }

	void MultiplyUInt64(uint64_t x, uint64_t y, uint64_t* prod_hi,
							uint64_t* prod_lo) {
	uint128_t prod = MultiplyUInt64(x, y);
	*prod_hi = static_cast<uint64_t>(prod >> 64);
	*prod_lo = static_cast<uint64_t>(prod);
	}

	uint128_t MultiplyUInt64(uint64_t x, uint64_t y) {
	return uint128_t(x) * uint128_t(y);
	}

	uint64_t MultiplyMod(uint64_t x, uint64_t y, uint64_t modulus) {

	if(modulus == 0) { std::cout << "modulus == 0" << std::endl; exit(0);}
	if(!(x < modulus)) { std::cout << "x " << x << " >= modulus " << modulus << std::endl; exit(0);}
	if(!(y < modulus)) { std::cout << "y " << y << " >= modulus " << modulus << std::endl; exit(0);}

	uint64_t prod_hi, prod_lo;
	MultiplyUInt64(x, y, &prod_hi, &prod_lo);

	return BarrettReduce128(prod_hi, prod_lo, modulus);
	}

	uint64_t BarrettReduce128(uint64_t input_hi, uint64_t input_lo,
									uint64_t modulus) {
	//   HEXL_CHECK(modulus != 0, "modulus == 0")
	if(modulus == 0) { std::cout << "modulus == 0" << std::endl; exit(0);}
	uint128_t n = (static_cast<uint128_t>(input_hi) << 64) |
					(static_cast<uint128_t>(input_lo));

	return n % modulus;
	// TODO(fboemer): actually use barrett reduction if performance-critical
	}

	uint64_t GenerateInsecureUniformRandomValue(uint64_t min_value,
													uint64_t max_value) {
	//   HEXL_CHECK(min_value < max_value, "min_value must be > max_value");

	static std::random_device rd;
	static std::mt19937 mersenne_engine(rd());
	std::uniform_int_distribution<uint64_t> distrib(min_value, max_value - 1);

	return distrib(mersenne_engine);
	}

	uint64_t Log2(uint64_t x) { return MSB(x); }

	uint64_t MSB(uint64_t input) {
	return static_cast<uint64_t>(std::log2l(input));
	}

	uint64_t GeneratePrimitiveRoot(uint64_t degree, uint64_t modulus) {
	// We need to divide modulus-1 by degree to get the size of the quotient group
	uint64_t size_entire_group = modulus - 1;

	// Compute size of quotient group
	uint64_t size_quotient_group = size_entire_group / degree;

	for (int trial = 0; trial < 200; ++trial) {
		uint64_t root = GenerateInsecureUniformRandomValue(0, modulus);
		root = PowMod(root, size_quotient_group, modulus);

		if (IsPrimitiveRoot(root, degree, modulus)) {
		return root;
		}
	}
	std::cout << "no primitive root found for degree "
	                        << degree << " modulus " << modulus << endl;
	return 0;
	}


	template <typename T>
	void print_vector(vector<T> vector_temp){
		cout << "{";
		for(auto item : vector_temp) cout << item << ",";
		cout << "\b}";
	}

	template <typename T>
	void print_array(T *array_temp,size_t len){
		cout << "{";
		for (size_t i = 0; i < len; i++)  cout << array_temp[i] << ",";
		cout << "\b}";
	}
};