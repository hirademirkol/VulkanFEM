#include "utils.hpp"

Vec3d operator*(Vec3i first, Vec3d second)
{
	Vec3d result;
	result.x = first.x * second.x;
	result.y = first.y * second.y;
	result.z = first.z * second.z;
	return result;
}

Vec3i operator+(Vec3i first, Vec3i second)
{
	Vec3i result;
	result.x = first.x + second.x;
	result.y = first.y + second.y;
	result.z = first.z + second.z;
	return result;
}

Vec3i operator-(Vec3i first, Vec3i second)
{
	Vec3i result;
	result.x = first.x - second.x;
	result.y = first.y - second.y;
	result.z = first.z - second.z;
	return result;
}

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