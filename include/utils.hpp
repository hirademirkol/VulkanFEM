#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <vector>

template <typename scalar>
struct Vec3
{
	scalar x;
	scalar y;
	scalar z;

	inline static scalar MAX;

	Vec3()
	{
		this->x = 0;
		this->y = 0;
		this->z = 0;
	}

	Vec3(scalar a)
	{
		this->x = a;
		this->y = a;
		this->z = a;
	}

	Vec3(scalar a, scalar b, scalar c)
	{
		this->x = a;
		this->y = b;
		this->z = c;
	}

	Vec3 operator+=(Vec3 &second)
	{
		Vec3 result;
		result->x += second.x;
		result->y += second.y;
		result->z += second.z;
		return result;
	}

	Vec3 operator-=(Vec3 &second)
	{
		Vec3 result;
		result->x -= second.x;
		result->y -= second.y;
		result->z -= second.z;
		return result;
	}

	Vec3 operator*=(Vec3 &second)
	{
		Vec3 result;
		result.x = this->x * second.x;
		result.y = this->y * second.y;
		result.z = this->z * second.z;
		return result;
	}

	bool operator==(const Vec3 &second) const
	{
		return this->x == second.x && this->y == second.y && this->z == second.z;
	}

	scalar multiplyComponents()
	{
		return this->x * this->y * this->z;
	}
};

typedef Vec3<int> Vec3i;
typedef Vec3<float> Vec3f;
typedef Vec3<double> Vec3d;

Vec3d operator*(Vec3i first, Vec3d second);
Vec3i operator+(Vec3i first, Vec3i second);
Vec3i operator-(Vec3i first, Vec3i second);
bool operator<(Vec3i first, Vec3i second);

#pragma region // TODO: Change these
template <typename T>
inline Vec3<T> Vectorize(const T &s, const T &base)
{
	return Vec3<T>(s % base, (s % (base * base)) / base, s / (base * base));
}
template <class T>
inline T Linearize(const Vec3<T> &v, const Vec3<T> &size) { return v.x + v.y * size.x + v.z * (size.x * size.y); }

#define FOR3(vi, vi0, vi1)											\
	for ((vi).z = (vi0).z; (vi).z < (vi1).z; (vi).z++)    \
		for ((vi).y = (vi0).y; (vi).y < (vi1).y; (vi).y++)\
			for ((vi).x = (vi0).x; (vi).x < (vi1).x; (vi).x++)

#define N_BITS_X 11
#define N_BITS_Y 11
#define N_BITS_Z 10
inline static int PackPosition(Vec3i vi)
{
	return vi.x | (vi.y << N_BITS_X) | (vi.z << (N_BITS_X + N_BITS_Y));
}

#pragma endregion

namespace std
{
	template <>
	struct hash<Vec3i>
	{
		std::size_t operator()(const Vec3i& k) const
		{
			return std::hash<int>()(Linearize(k, Vec3i(k.MAX+1)));
		}
	};
} // namespace std

#endif // __UTILS_HPP__