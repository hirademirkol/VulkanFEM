#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <vector>

// Utility Vector (x,y,z) struct for basic coordinate operations
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

// Linearization of 3-components into 1 for node ordering
template <class T>
inline T Linearize(const Vec3<T> &v, const Vec3<T> &size) { return v.x + v.y * size.x + v.z * (size.x * size.y); }

#define FOR3(vi, vi0, vi1)								  \
	for ((vi).z = (vi0).z; (vi).z < (vi1).z; (vi).z++)    \
		for ((vi).y = (vi0).y; (vi).y < (vi1).y; (vi).y++)\
			for ((vi).x = (vi0).x; (vi).x < (vi1).x; (vi).x++)

namespace std
{
	// Hashing operator for finding unique coordinates
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