/*
Copyright (c) 2023, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#ifndef POISSON_RECON_INCLUDED
#define POISSON_RECON_INCLUDED

#include "FEMTree.h"
#include "VertexFactory.h"

#undef USE_DOUBLE								// If enabled, double-precesion is used

#define DATA_DEGREE 0							// The order of the B-Spline used to splat in data for color interpolation
#define WEIGHT_DEGREE 2							// The order of the B-Spline used to splat in the weights for density estimation
#define NORMAL_DEGREE 2							// The order of the B-Spline used to splat in the normals for constructing the Laplacian constraints
#define DEFAULT_FEM_DEGREE 1					// The default finite-element degree
#define DEFAULT_FEM_BOUNDARY BOUNDARY_NEUMANN	// The default finite-element boundary type {BOUNDARY_FREE, BOUNDARY_DIRICHLET, BOUNDARY_NEUMANN}
#define DEFAULT_DIMENSION 3						// The dimension of the system

const float DefaultPointWeightMultiplier = 2.f;

template< unsigned int Dim , typename Real >
struct ConstraintDual
{
	Real target , weight;
	ConstraintDual( void ) : target(0) , weight(0) {}
	ConstraintDual( Real t , Real w ) : target(t) , weight(w) {}
	CumulativeDerivativeValues< Real , Dim , 0 > operator()( const Point< Real , Dim >& p ) const { return CumulativeDerivativeValues< Real , Dim , 0 >( target*weight ); };
};

template< unsigned int Dim , typename Real >
struct SystemDual
{
	Real weight;
	SystemDual( void ) : weight(0) {}
	SystemDual( Real w ) : weight(w) {}
	CumulativeDerivativeValues< Real , Dim , 0 > operator()( const Point< Real , Dim >& p , const CumulativeDerivativeValues< Real , Dim , 0 >& dValues ) const { return dValues * weight; };
	CumulativeDerivativeValues< double , Dim , 0 > operator()( const Point< Real , Dim >& p , const CumulativeDerivativeValues< double , Dim , 0 >& dValues ) const { return dValues * weight; };
};

template< unsigned int Dim >
struct SystemDual< Dim , double >
{
	typedef double Real;
	Real weight;
	SystemDual( void ) : weight(0) {}
	SystemDual( Real w ) : weight(w) {}
	CumulativeDerivativeValues< Real , Dim , 0 > operator()( const Point< Real , Dim >& p , const CumulativeDerivativeValues< Real , Dim , 0 >& dValues ) const { return dValues * weight; };
};

template< typename Real , unsigned int Dim , typename AuxData >
using InputOrientedPointStreamInfo = typename FEMTreeInitializer< Dim , Real >::template InputPointStream< VectorTypeUnion< Real , typename VertexFactory::NormalFactory< Real , Dim >::VertexType , AuxData > >;

template< typename Real , unsigned int Dim , typename AuxData >
using InputOrientedPointStream = typename InputOrientedPointStreamInfo< Real , Dim , AuxData >::StreamType;

template< typename Real , unsigned int Dim , typename AuxData >
using OutputOrientedPointStreamInfo = typename FEMTreeInitializer< Dim , Real >::template OutputPointStream< VectorTypeUnion< Real , typename VertexFactory::NormalFactory< Real , Dim >::VertexType , AuxData > >;

template< typename Real , unsigned int Dim , typename AuxData >
using OutputOrientedPointStream = typename OutputOrientedPointStreamInfo< Real , Dim , AuxData >::StreamType;

template< class Real , unsigned int Dim >
XForm< Real , Dim+1 > GetBoundingBoxXForm( Point< Real , Dim > min , Point< Real , Dim > max , Real scaleFactor )
{
	Point< Real , Dim > center = ( max + min ) / 2;
	Real scale = max[0] - min[0];
	for( int d=1 ; d<Dim ; d++ ) scale = std::max< Real >( scale , max[d]-min[d] );
	scale *= scaleFactor;
	for( int i=0 ; i<Dim ; i++ ) center[i] -= scale/2;
	XForm< Real , Dim+1 > tXForm = XForm< Real , Dim+1 >::Identity() , sXForm = XForm< Real , Dim+1 >::Identity();
	for( int i=0 ; i<Dim ; i++ ) sXForm(i,i) = (Real)(1./scale ) , tXForm(Dim,i) = -center[i];
	unsigned int maxDim = 0;
	for( int i=1 ; i<Dim ; i++ ) if( (max[i]-min[i])>(max[maxDim]-min[maxDim]) ) maxDim = i;
	XForm< Real , Dim+1 > rXForm;
	for( int i=0 ; i<Dim ; i++ ) rXForm((maxDim+i)%Dim,(Dim-1+i)%Dim) = 1;
	rXForm(Dim,Dim) = 1;
	return rXForm * sXForm * tXForm;
}

template< class Real , unsigned int Dim , typename AuxData >
XForm< Real , Dim+1 > GetPointXForm( InputOrientedPointStream< Real , Dim , AuxData > &stream , typename InputOrientedPointStreamInfo< Real , Dim , AuxData >::DataType d , Real scaleFactor )
{
	Point< Real , Dim > min , max;
	InputOrientedPointStreamInfo< Real , Dim , AuxData >::BoundingBox( stream , d , min , max );
	return GetBoundingBoxXForm( min , max , scaleFactor );
}

#endif // POISSON_RECON_INCLUDED