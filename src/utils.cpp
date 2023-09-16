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

bool operator<(Vec3i first, Vec3i second)
{
	return Linearize(first, Vec3i(first.MAX+1)) < Linearize(second, Vec3i(second.MAX+1));
}