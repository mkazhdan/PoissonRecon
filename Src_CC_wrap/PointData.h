//##########################################################################
//#                                                                        #
//#               CLOUDCOMPARE WRAPPER: PoissonReconLib                    #
//#                                                                        #
//#  This program is free software; you can redistribute it and/or modify  #
//#  it under the terms of the GNU General Public License as published by  #
//#  the Free Software Foundation; version 2 or later of the License.      #
//#                                                                        #
//#  This program is distributed in the hope that it will be useful,       #
//#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
//#  GNU General Public License for more details.                          #
//#                                                                        #
//#               COPYRIGHT: Daniel Girardeau-Montaut                      #
//#                                                                        #
//##########################################################################

#ifndef POINTDATA_H
#define POINTDATA_H

template <typename Real>
class PointData {
public:
	PointData() : normal{ 0, 0, 0 }, color{ 0, 0, 0 } {}
	PointData(const Real _normal[3], const Real _color[3], Real scale = 1.0)
	{
		normal[0] = scale * _normal[0];
		normal[1] = scale * _normal[1];
		normal[2] = scale * _normal[2];
		color[0] = scale * _color[0];
		color[1] = scale * _color[1];
		color[2] = scale * _color[2];
	}

	PointData operator * (Real s) const
	{
		return PointData(normal, color, s);
	}

	PointData operator / (Real s) const
	{
		return PointData(normal, color, 1 / s);
	}

	PointData& operator += (const PointData& d)
	{
		normal[0] += d.normal[0];
		normal[1] += d.normal[1];
		normal[2] += d.normal[2];
		color[0] += d.color[0];
		color[1] += d.color[1];
		color[2] += d.color[2];
		return *this;
	}
	PointData& operator *= (Real s)
	{
		normal[0] *= s;
		normal[1] *= s;
		normal[2] *= s;
		color[0] *= s;
		color[1] *= s;
		color[2] *= s;
		return *this;
	}

public:
	Real normal[3];
	Real color[3];
};

extern template class PointData<float>;
extern template class PointData<double>;

#endif // POINTDATA_H
