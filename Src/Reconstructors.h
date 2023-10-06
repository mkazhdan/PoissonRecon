/*
Copyright (c) 2022, Michael Kazhdan and Matthew Bolitho
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

#ifndef RECONSTRUCTORS_INCLUDED
#define RECONSTRUCTORS_INCLUDED

#include "MyMiscellany.h"
#include "DataStream.imp.h"
#include "FEMTree.h"

namespace Reconstructor
{
	static const unsigned int DataDegree = 0;	// The order of the B-Spline used to splat in data for auxiliary data interpolation
	static const unsigned int WeightDegree = 2;	// The order of the B-Spline used to splat in the weights for density estimation

	// For clarity, to distinguies betwen the case that a Point is referencing a position and a normal
	template< typename Real , unsigned int Dim > using Position = Point< Real , Dim >;
	template< typename Real , unsigned int Dim > using Normal = Point< Real , Dim >;
	template< typename Real , unsigned int Dim > using Gradient = Point< Real , Dim >;

#include "Reconstructors.streams.h"

	// Declare a type for storing the solution information
	template< typename Real , unsigned int Dim , unsigned int FEMSig , typename ... AuxData > struct ReconstructionInfo;

	// Specialized solution information without auxiliary data
	template< typename Real , unsigned int Dim , unsigned int FEMSig >
	struct ReconstructionInfo< Real , Dim , FEMSig >
	{
		// The signature pack
		typedef IsotropicUIntPack< Dim , FEMSig > Sigs;

		// The type representing the point sampling density
		typedef typename FEMTree< Dim , Real >::template DensityEstimator< Reconstructor::WeightDegree > DensityEstimator;

		// The constructor
		ReconstructionInfo( void ) : density(NULL) , isoValue(0) , tree(MEMORY_ALLOCATOR_BLOCK_SIZE) , unitCubeToModel( XForm< Real , Dim+1 >::Identity() ){}

		// The desctructor
		~ReconstructionInfo( void ){ delete density ; density = NULL; }

		// The transformation taking points in the unit cube back to world coordinates
		XForm< Real , Dim+1 > unitCubeToModel;

		// The octree adapted to the points
		FEMTree< Dim , Real > tree;

		// The solution coefficients
		DenseNodeData< Real , Sigs > solution;

		// The average value at the sample positions
		Real isoValue;

		// The density estimator computed from the samples
		DensityEstimator *density;
	};

	// Specialized solution information with auxiliary data
	template< typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData >
	struct ReconstructionInfo< Real , Dim , FEMSig , AuxData > : public ReconstructionInfo< Real , Dim , FEMSig >
	{
		typedef IsotropicUIntPack< Dim , FEMSig > Sigs;

		// The signature of the finite-element used for data extrapolation
		static const unsigned int DataSig = FEMDegreeAndBType< Reconstructor::DataDegree , BOUNDARY_FREE >::Signature;

		// The constructor
		ReconstructionInfo( AuxData zeroAuxData ) : auxData(NULL) , zeroAuxData(zeroAuxData) {}

		// The desctructor
		~ReconstructionInfo( void ){ delete auxData ; auxData = NULL; }

		// The auxiliary information stored with the oriented vertices
		SparseNodeData< ProjectiveData< AuxData , Real > , IsotropicUIntPack< Dim , DataSig > > *auxData;

		// An instance of "zero" AuxData
		AuxData zeroAuxData;

		// a method for rescaling the contents of the auxiliary data to give interpolation-preference to finer levels
		void weightAuxDataByDepth( Real perLevelScaleFactor )
		{
			auto nodeFunctor = [&]( const RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type > *n )
			{
				ProjectiveData< AuxData , Real >* clr = (*auxData)( n );
				if( clr ) (*clr) *= (Real)pow( (Real)perLevelScaleFactor , ReconstructionInfo< Real , Dim , FEMSig>::tree.depth( n ) );
			};
			ReconstructionInfo< Real , Dim , FEMSig>::tree.tree().processNodes( nodeFunctor );
		}
	};

	struct MeshExtractionParameters
	{
		bool linearFit;
		bool outputGradients;
		bool forceManifold;
		bool polygonMesh;
		bool verbose;
		MeshExtractionParameters( void ) : linearFit(false) , outputGradients(false) , forceManifold(true) , polygonMesh(false) , verbose(false) {}
	};

	template< bool HasAuxData , typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData , typename OutputVertexStream , typename ReconstructionInfoType , unsigned int ... FEMSigs >
	void _ExtractMesh( UIntPack< FEMSigs ... >  , const ReconstructionInfoType &sInfo , OutputVertexStream &vertexStream , OutputDataStream< std::vector< node_index_type > > &polygonStream , MeshExtractionParameters params );

	template< typename Real , unsigned int Dim , unsigned int FEMSig >
	void ExtractMesh( const ReconstructionInfo< Real , Dim , FEMSig > &sInfo , OutputVertexStream< Real , Dim > &vertexStream , OutputDataStream< std::vector< node_index_type > > &polygonStream , MeshExtractionParameters params )
	{
		typedef unsigned char AuxData;
		_ExtractMesh< false , Real , Dim , FEMSig , AuxData , OutputVertexStream< Real , Dim > , ReconstructionInfo< Real , Dim , FEMSig > >( IsotropicUIntPack< Dim , FEMSig >() , sInfo , vertexStream , polygonStream , params );
	}

	template< typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData >
	void ExtractMesh( const ReconstructionInfo< Real , Dim , FEMSig , AuxData > &sInfo , OutputVertexWithDataStream< Real , Dim , AuxData > &vertexStream , OutputDataStream< std::vector< node_index_type > > &polygonStream , MeshExtractionParameters params )
	{
		_ExtractMesh< true , Real , Dim , FEMSig , AuxData , OutputVertexWithDataStream< Real , Dim , AuxData > , ReconstructionInfo< Real , Dim , FEMSig , AuxData > >( IsotropicUIntPack< Dim , FEMSig >() , sInfo , vertexStream , polygonStream , params );
	}

	namespace Poisson
	{
		static const unsigned int NormalDegree = 2;							// The order of the B-Spline used to splat in the normals for constructing the Laplacian constraints
		static const unsigned int DefaultFEMDegree = 1;						// The default finite-element degree (has to be at least 1)
		static const BoundaryType DefaultFEMBoundary = BOUNDARY_NEUMANN;	// The default finite-element boundary type {BOUNDARY_FREE, BOUNDARY_DIRICHLET, BOUNDARY_NEUMANN}
		static const float WeightMultiplier = 2.f;							// The default degree-to-point-weight scaling

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

		template< typename Real >
		struct SolutionParameters
		{
			bool verbose;
			bool dirichletErode;
			bool outputDensity;
			bool exactInterpolation;
			bool showResidual;
			Real scale;
			Real confidence;
			Real confidenceBias;
			Real lowDepthCutOff;
			Real width;
			Real pointWeight;
			Real samplesPerNode;
			Real cgSolverAccuracy;
#ifdef SOFT_DIRICHLET
			Real dirichletWeight;
#endif // SOFT_DIRICHLET
			unsigned int depth;
			unsigned int solveDepth;
			unsigned int baseDepth;
			unsigned int fullDepth;
			unsigned int kernelDepth;
			unsigned int envelopeDepth;
			unsigned int baseVCycles;
			unsigned int iters;

			SolutionParameters( void ) :
				verbose(false) , dirichletErode(false) , outputDensity(false) , exactInterpolation(false) , showResidual(false) ,
				scale((Real)1.1) , confidence((Real)0.) , confidenceBias((Real)0.) , lowDepthCutOff((Real)0.) , width((Real)0.) ,
				pointWeight((Real)0.) , samplesPerNode((Real)1.5) , cgSolverAccuracy((Real)1e-3 ) ,
#ifdef SOFT_DIRICHLET
				dirichletWeight((Real)0.) ,
#endif // SOFT_DIRICHLET
				depth((unsigned int)8) , solveDepth((unsigned int)-1) , baseDepth((unsigned int)-1) , fullDepth((unsigned int)5) , kernelDepth((unsigned int)-1) ,
				envelopeDepth((unsigned int)-1) , baseVCycles((unsigned int)1) , iters((unsigned int)8)
			{}
		};

		template< typename Real , unsigned int Dim >
		struct EnvelopeMesh
		{
			std::vector< Point< Real , Dim > > vertices;
			std::vector< SimplexIndex< Dim-1 , node_index_type > > simplices;
		};

		template< bool HasAuxData , typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData , typename InputSampleStreamType , unsigned int ... FEMSigs >
		static typename std::conditional< HasAuxData , ReconstructionInfo< Real , Dim , FEMSig , AuxData > , ReconstructionInfo< Real , Dim , FEMSig > >::type *_Solve( UIntPack< FEMSigs... > , InputSampleStreamType &pointStream , SolutionParameters< Real > params , const EnvelopeMesh< Real , Dim > *envelopeMesh );

#ifdef DE_VIRTUALIZE_INPUT
		template< typename Real , unsigned int Dim , unsigned int FEMSig , typename InputSampleStreamType >
		ReconstructionInfo< Real , Dim , FEMSig > *Solve( InputSampleStreamType &pointStream , SolutionParameters< Real > params , const EnvelopeMesh< Real , Dim > *envelopeMesh=NULL )
		{
			static_assert( std::is_base_of< InputSampleStream< Real , Dim > , InputSampleStreamType >::value , "[ERROR] Unexpected sample stream type" );
			typedef unsigned char AuxData;
			return _Solve< false , Real , Dim , FEMSig , AuxData , InputSampleStreamType >( IsotropicUIntPack< Dim , FEMSig >() , pointStream , params , envelopeMesh );
		}

		template< typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData , typename InputSampleStreamType >
		ReconstructionInfo< Real , Dim , FEMSig , AuxData > *Solve( InputSampleStreamType &pointStream , SolutionParameters< Real > params , const EnvelopeMesh< Real , Dim > *envelopeMesh=NULL )
		{
			static_assert( std::is_base_of< InputSampleWithDataStream< Real , Dim , AuxData > , InputSampleStreamType >::value , "[ERROR] Unexpected sample stream type" );
			return _Solve< true , Real , Dim , FEMSig , AuxData , InputSampleStreamType >( IsotropicUIntPack< Dim , FEMSig >() , pointStream , params , envelopeMesh );
		}
#else // !DE_VIRTUALIZE_INPUT
		template< typename Real , unsigned int Dim , unsigned int FEMSig >
		ReconstructionInfo< Real , Dim , FEMSig > *Solve( InputSampleStream< Real , Dim > &pointStream , SolutionParameters< Real > params , const EnvelopeMesh< Real , Dim > *envelopeMesh=NULL )
		{
			typedef unsigned char AuxData;
			return _Solve< false , Real , Dim , FEMSig , AuxData , InputSampleStream< Real , Dim > >( IsotropicUIntPack< Dim , FEMSig >() , pointStream , params , envelopeMesh );
		}

		template< typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData >
		ReconstructionInfo< Real , Dim , FEMSig , AuxData > *Solve( InputSampleWithDataStream< Real , Dim , AuxData > &pointStream , SolutionParameters< Real > params , const EnvelopeMesh< Real , Dim > *envelopeMesh=NULL )
		{
			return _Solve< true , Real , Dim , FEMSig , AuxData , InputSampleWithDataStream< Real , Dim , AuxData > >( IsotropicUIntPack< Dim , FEMSig >() , pointStream , params , envelopeMesh );
		}
#endif // DE_VIRTUALIZE_INPUT
	}

	namespace SSD
	{
		static const unsigned int NormalDegree = 2;								// The order of the B-Spline used to splat in the normals for constructing the Laplacian constraints
		static const unsigned int DefaultFEMDegree = 2;							// The default finite-element degree (has to be at least 2)
		static const BoundaryType DefaultFEMBoundary = BOUNDARY_NEUMANN;		// The default finite-element boundary type {BOUNDARY_FREE, BOUNDARY_DIRICHLET, BOUNDARY_NEUMANN}
		static const double WeightMultipliers[] = { 5e+1f , 5e-4f , 1e-5f };	// The default weights for balancing the value, gradient, and laplacian energy terms

		template< unsigned int Dim , typename ... > struct ConstraintDual;
		template< unsigned int Dim , typename ... > struct SystemDual;

		template< unsigned int Dim , typename Real >
		struct ConstraintDual< Dim , Real >
		{
			Real target , vWeight , gWeight;
			ConstraintDual( Real t , Real v , Real g ) : target(t) , vWeight(v) , gWeight(g) { }
			CumulativeDerivativeValues< Real , Dim , 1 > operator()( const Point< Real , Dim >& p , const Point< Real , Dim > &n ) const 
			{
				CumulativeDerivativeValues< Real , Dim , 1 > cdv;
				cdv[0] = target*vWeight;
				for( int d=0 ; d<Dim ; d++ ) cdv[1+d] = -n[d]*gWeight;
				return cdv;
			}
		};

		template< unsigned int Dim , typename Real , typename AuxData >
		struct ConstraintDual< Dim , Real , AuxData >
		{
			Real target , vWeight , gWeight;
			ConstraintDual( Real t , Real v , Real g ) : target(t) , vWeight(v) , gWeight(g) { }
			CumulativeDerivativeValues< Real , Dim , 1 > operator()( const Point< Real , Dim >& p , const VectorTypeUnion< Real , Point< Real , Dim > , AuxData > &normalAndAuxData ) const 
			{
				Point< Real , Dim > n = normalAndAuxData.template get<0>();
				CumulativeDerivativeValues< Real , Dim , 1 > cdv;
				cdv[0] = target*vWeight;
				for( int d=0 ; d<Dim ; d++ ) cdv[1+d] = -n[d]*gWeight;
				return cdv;
			}
		};

		template< unsigned int Dim , typename Real >
		struct SystemDual< Dim , Real >
		{
			CumulativeDerivativeValues< Real , Dim , 1 > weight;
			SystemDual( Real v , Real g )
			{
				weight[0] = v;
				for( int d=0 ; d<Dim ; d++ ) weight[d+1] = g;
			}
			CumulativeDerivativeValues< Real , Dim , 1 > operator()( Point< Real , Dim > p , const Point< Real , Dim > & , const CumulativeDerivativeValues< Real , Dim , 1 >& dValues ) const
			{
				return dValues * weight;
			}
			CumulativeDerivativeValues< double , Dim , 1 > operator()( Point< Real , Dim > p , const Point< Real , Dim > & , const CumulativeDerivativeValues< double , Dim , 1 >& dValues ) const
			{
				return dValues * weight;
			};
		};

		template< unsigned int Dim >
		struct SystemDual< Dim , double >
		{
			typedef double Real;
			CumulativeDerivativeValues< Real , Dim , 1 > weight;
			SystemDual( Real v , Real g ) : weight( v , g , g , g ) { }
			CumulativeDerivativeValues< Real , Dim , 1 > operator()( Point< Real , Dim > p , const Point< Real , Dim > & , const CumulativeDerivativeValues< Real , Dim , 1 >& dValues ) const
			{
				return dValues * weight;
			}
		};

		template< unsigned int Dim , typename Real , typename AuxData >
		struct SystemDual< Dim , Real , AuxData >
		{
			CumulativeDerivativeValues< Real , Dim , 1 > weight;
			SystemDual( Real v , Real g )
			{
				weight[0] = v;
				for( int d=0 ; d<Dim ; d++ ) weight[d+1] = g;
			}
			CumulativeDerivativeValues< Real , Dim , 1 > operator()( Point< Real , Dim > p , const VectorTypeUnion< Real , Point< Real , Dim > , AuxData > & , const CumulativeDerivativeValues< Real , Dim , 1 >& dValues ) const
			{
				return dValues * weight;
			}
			CumulativeDerivativeValues< double , Dim , 1 > operator()( Point< Real , Dim > p , const VectorTypeUnion< Real , Point< Real , Dim > , AuxData > & , const CumulativeDerivativeValues< double , Dim , 1 >& dValues ) const
			{
				return dValues * weight;
			};
		};

		template< unsigned int Dim , class AuxData >
		struct SystemDual< Dim , double , AuxData >
		{
			typedef double Real;
			CumulativeDerivativeValues< Real , Dim , 1 > weight;
			SystemDual( Real v , Real g ) : weight( v , g , g , g ) { }
			CumulativeDerivativeValues< Real , Dim , 1 > operator()( Point< Real , Dim > p , const VectorTypeUnion< Real , Point< Real , Dim > , AuxData > & , const CumulativeDerivativeValues< Real , Dim , 1 >& dValues ) const
			{
				return dValues * weight;
			}
		};

		template< typename Real >
		struct SolutionParameters
		{
			bool verbose;
			bool outputDensity;
			bool exactInterpolation;
			bool showResidual;
			Real scale;
			Real confidence;
			Real confidenceBias;
			Real lowDepthCutOff;
			Real width;
			Real pointWeight;
			Real gradientWeight;
			Real biLapWeight;
			Real samplesPerNode;
			Real cgSolverAccuracy;
			unsigned int depth;
			unsigned int solveDepth;
			unsigned int baseDepth;
			unsigned int fullDepth;
			unsigned int kernelDepth;
			unsigned int baseVCycles;
			unsigned int iters;

			SolutionParameters( void ) :
				verbose(false) , outputDensity(false) , exactInterpolation(false) , showResidual(false) ,
				scale((Real)1.1) , confidence((Real)0.) , confidenceBias((Real)0.) , lowDepthCutOff((Real)0.) , width((Real)0.) ,
				pointWeight((Real)WeightMultipliers[0]) , gradientWeight((Real)WeightMultipliers[1]) , biLapWeight((Real)WeightMultipliers[2]) , samplesPerNode((Real)1.5) , cgSolverAccuracy((Real)1e-3 ) ,
				depth((unsigned int)8) , solveDepth((unsigned int)-1) , baseDepth((unsigned int)-1) , fullDepth((unsigned int)5) , kernelDepth((unsigned int)-1) ,
				baseVCycles((unsigned int)1) , iters((unsigned int)8)
			{}

		};

		template< bool HasAuxData , typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData , typename InputSampleStreamType , unsigned int ... FEMSigs >
		typename std::conditional< HasAuxData , ReconstructionInfo< Real , Dim , FEMSig , AuxData > , ReconstructionInfo< Real , Dim , FEMSig > >::type *_Solve( UIntPack< FEMSigs ... > , InputSampleStreamType &pointStream , SolutionParameters< Real > params );

#ifdef DE_VIRTUALIZE_INPUT
		template< typename Real , unsigned int Dim , unsigned int FEMSig , typename InputSampleStreamType >
		ReconstructionInfo< Real , Dim , FEMSig > *Solve( InputSampleStreamType &pointStream , SolutionParameters< Real > params )
		{
			static_assert( std::is_base_of< InputSampleStream< Real , Dim > , InputSampleStreamType >::value , "[ERROR] Unexpected sample stream type" );
			typedef unsigned char AuxData;
			return _Solve< false , Real , Dim , FEMSig , AuxData , InputSampleStreamType >( IsotropicUIntPack< Dim , FEMSig >() , pointStream , params );
		}

		template< typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData , typename InputSampleStreamType >
		ReconstructionInfo< Real , Dim , FEMSig , AuxData > *Solve( InputSampleStreamType &pointStream , SolutionParameters< Real > params )
		{
			static_assert( std::is_base_of< InputSampleWithDataStream< Real , Dim , AuxData > , InputSampleStreamType >::value , "[ERROR] Unexpected sample stream type" );
			return _Solve< true , Real , Dim , FEMSig , AuxData , InputSampleStreamType >( IsotropicUIntPack< Dim , FEMSig >() , pointStream , params );
		}
#else // !DE_VIRTUALIZE_INPUT
		template< typename Real , unsigned int Dim , unsigned int FEMSig >
		ReconstructionInfo< Real , Dim , FEMSig > *Solve( InputSampleStream< Real , Dim > &pointStream , SolutionParameters< Real > params )
		{
			typedef unsigned char AuxData;
			return _Solve< false , Real , Dim , FEMSig , AuxData , InputSampleStream< Real , Dim > >( IsotropicUIntPack< Dim , FEMSig >() , pointStream , params );
		}

		template< typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData >
		ReconstructionInfo< Real , Dim , FEMSig , AuxData > *Solve( InputSampleWithDataStream< Real , Dim , AuxData > &pointStream , SolutionParameters< Real > params )
		{
			return _Solve< true , Real , Dim , FEMSig , AuxData , InputSampleWithDataStream< Real , Dim , AuxData > >( IsotropicUIntPack< Dim , FEMSig >() , pointStream , params );
		}
#endif // DE_VIRTUALIZE_INPUT
	}

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

#ifdef DE_VIRTUALIZE_INPUT
	template< class Real , unsigned int Dim , typename SampleStream >
	void SetBoundingBox( SampleStream &stream , Point< Real , Dim >& min , Point< Real , Dim >& max )
	{
		using Sample = Point< Real , Dim >;
		static_assert( std::is_base_of< InputDataStream< Sample > , SampleStream >::value , "[ERROR] Unexpected sample stream type" );
		Sample s;
		for( unsigned int d=0 ; d<Dim ; d++ ) min[d] = std::numeric_limits< Real >::infinity() , max[d] = -std::numeric_limits< Real >::infinity();
		while( stream.read( s ) ) for( unsigned int d=0 ; d<Dim ; d++ ) min[d] = std::min< Real >( min[d] , s[d] ) , max[d] = std::max< Real >( max[d] , s[d] );
		stream.reset();
	}

	template< class Real , unsigned int Dim , typename AuxData , typename SampleStream >
	void SetBoundingBox( SampleStream &stream , AuxData d , Point< Real , Dim >& min , Point< Real , Dim >& max )
	{
		using Sample = VectorTypeUnion< Real , Point< Real , Dim > , AuxData >;
		static_assert( std::is_base_of< InputDataStream< Sample > , SampleStream >::value , "[ERROR] Unexpected sample stream type" );
		Sample s( Point< Real , Dim >() , d );
		for( unsigned int d=0 ; d<Dim ; d++ ) min[d] = std::numeric_limits< Real >::infinity() , max[d] = -std::numeric_limits< Real >::infinity();
		while( stream.read( s ) ) for( unsigned int d=0 ; d<Dim ; d++ ) min[d] = std::min< Real >( min[d] , s.template get<0>()[d] ) , max[d] = std::max< Real >( max[d] , s.template get<0>()[d] );
		stream.reset();
	}

	template< class Real , unsigned int Dim , typename SampleStream >
	XForm< Real , Dim+1 > GetPointXForm( SampleStream &stream , Real scaleFactor )
	{
		using Sample = Point< Real , Dim >;
		static_assert( std::is_base_of< InputDataStream< Sample > , SampleStream >::value , "[ERROR] Unexpected sample stream type" );
		Point< Real , Dim > min , max;
		SetBoundingBox< Real , Dim , SampleStream >( stream , min , max );
		return GetBoundingBoxXForm( min , max , scaleFactor );
	}

	template< class Real , unsigned int Dim , typename AuxData , typename SampleStream >
	XForm< Real , Dim+1 > GetPointXForm( SampleStream &stream , AuxData d , Real scaleFactor )
	{
		using Sample = VectorTypeUnion< Real , Point< Real , Dim > , AuxData >;
		static_assert( std::is_base_of< InputDataStream< Sample > , SampleStream >::value , "[ERROR] Unexpected sample stream type" );
		Point< Real , Dim > min , max;
		SetBoundingBox< Real , Dim , AuxData , SampleStream >( stream , d , min , max );
		return GetBoundingBoxXForm( min , max , scaleFactor );
	}
#else // !DE_VIRTUALIZE_INPUT
	template< class Real , unsigned int Dim >
	void SetBoundingBox( InputDataStream< Point< Real , Dim > > &stream , Point< Real , Dim >& min , Point< Real , Dim >& max )
	{
		Point< Real , Dim > p;
		for( unsigned int d=0 ; d<Dim ; d++ ) min[d] = std::numeric_limits< Real >::infinity() , max[d] = -std::numeric_limits< Real >::infinity();
		while( stream.read( p ) ) for( unsigned int d=0 ; d<Dim ; d++ ) min[d] = std::min< Real >( min[d] , p[d] ) , max[d] = std::max< Real >( max[d] , p[d] );
		stream.reset();
	}

	template< class Real , unsigned int Dim , typename AuxData >
	void SetBoundingBox( InputDataStream< VectorTypeUnion< Real , Point< Real , Dim > , AuxData > > &stream , AuxData d , Point< Real , Dim >& min , Point< Real , Dim >& max )
	{
		VectorTypeUnion< Real , Point< Real , Dim > , AuxData > p( Point< Real , Dim >() , d );
		for( unsigned int d=0 ; d<Dim ; d++ ) min[d] = std::numeric_limits< Real >::infinity() , max[d] = -std::numeric_limits< Real >::infinity();
		while( stream.read( p ) ) for( unsigned int d=0 ; d<Dim ; d++ ) min[d] = std::min< Real >( min[d] , p.template get<0>()[d] ) , max[d] = std::max< Real >( max[d] , p.template get<0>()[d] );
		stream.reset();
	}

	template< class Real , unsigned int Dim >
	XForm< Real , Dim+1 > GetPointXForm( InputDataStream< Point< Real , Dim > > &stream , Real scaleFactor )
	{
		Point< Real , Dim > min , max;
		SetBoundingBox( stream , min , max );
		return GetBoundingBoxXForm( min , max , scaleFactor );
	}

	template< class Real , unsigned int Dim , typename AuxData >
	XForm< Real , Dim+1 > GetPointXForm( InputDataStream< VectorTypeUnion< Real , Point< Real , Dim > , AuxData > > &stream , AuxData d , Real scaleFactor )
	{
		Point< Real , Dim > min , max;
		SetBoundingBox( stream , d , min , max );
		return GetBoundingBoxXForm( min , max , scaleFactor );
	}
#endif // DE_VIRTUALIZE_INPUT

	template< bool HasAuxData , typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData , typename OutputVertexStream , typename ReconstructionInfoType , unsigned int ... FEMSigs >
	void _ExtractMesh( UIntPack< FEMSigs ... > , const ReconstructionInfoType &sInfo , OutputVertexStream &vertexStream , OutputDataStream< std::vector< node_index_type > > &polygonStream , MeshExtractionParameters params )
	{
		typedef UIntPack< FEMSigs ... > Sigs;
		static_assert( std::is_same< IsotropicUIntPack< Dim , FEMSig > , UIntPack< FEMSigs ... > >::value , "[ERROR] Signatures don't match" );
		static const unsigned int DataSig = FEMDegreeAndBType< Reconstructor::DataDegree , BOUNDARY_FREE >::Signature;
		typedef typename ReconstructionInfoType::DensityEstimator DensityEstimator;

		if constexpr( Dim!=3 )
		{
			WARN( "Extraction only supported for dimension 3" );
			return;
		}
		else
		{
			Profiler profiler(20);

			std::string statsString;

			if constexpr( HasAuxData )
			{
				typename LevelSetExtractor< Real , Dim , AuxData >::Stats stats;
				TransformedOutputVertexWithDataStream< Real , Dim , AuxData > _vertexStream( sInfo.unitCubeToModel , vertexStream );
				stats = LevelSetExtractor< Real , Dim , AuxData >::Extract( Sigs() , UIntPack< Reconstructor::WeightDegree >() , UIntPack< DataSig >() , sInfo.tree , sInfo.density , sInfo.auxData , sInfo.solution , sInfo.isoValue , _vertexStream , polygonStream , sInfo.zeroAuxData , !params.linearFit , params.outputGradients , params.forceManifold , params.polygonMesh , false );
				statsString = stats.toString();
			}
			else
			{
				typename LevelSetExtractor< Real , Dim >::Stats stats;
				TransformedOutputVertexStream< Real , Dim > _vertexStream( sInfo.unitCubeToModel , vertexStream );
				stats = LevelSetExtractor< Real , Dim >::Extract( Sigs() , UIntPack< Reconstructor::WeightDegree >() , sInfo.tree , sInfo.density , sInfo.solution , sInfo.isoValue , _vertexStream , polygonStream , !params.linearFit , params.outputGradients , params.forceManifold , params.polygonMesh , false );
				statsString = stats.toString();
			}
			if( params.verbose )
			{
				std::cout << "Vertices / Polygons: " << vertexStream.size() << " / " << polygonStream.size() << std::endl;
				std::cout << statsString << std::endl;
				if( params.polygonMesh ) std::cout << "#         Got polygons: " << profiler << std::endl;
				else                     std::cout << "#        Got triangles: " << profiler << std::endl;
			}
		}
	}

	template< bool HasAuxData , typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData , typename InputSampleStreamType , unsigned int ... FEMSigs >
	static typename std::conditional< HasAuxData , ReconstructionInfo< Real , Dim , FEMSig , AuxData > , ReconstructionInfo< Real , Dim , FEMSig > >::type *Poisson::_Solve( UIntPack< FEMSigs ... > , InputSampleStreamType &pointStream , SolutionParameters< Real > params , const EnvelopeMesh< Real , Dim > *envelopeMesh )
	{
		static_assert( std::is_same< IsotropicUIntPack< Dim , FEMSig > , UIntPack< FEMSigs... > >::value , "[ERROR] Signatures don't match" );

		// The signature for the finite-elements representing the auxiliary data (if it's there)
		static const unsigned int DataSig = FEMDegreeAndBType< Reconstructor::DataDegree , BOUNDARY_FREE >::Signature;

		///////////////
		// Types --> //
		// The packed finite elements signature
		typedef UIntPack< FEMSigs ... > Sigs;

		// The degrees of the finite elements across the different axes
		typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > Degrees;

		// The signature describing the normals elements
		typedef UIntPack< FEMDegreeAndBType< Poisson::NormalDegree , DerivativeBoundary< FEMSignature< FEMSigs >::BType , 1 >::BType >::Signature ... > NormalSigs;

		// Type for tracking sample interpolation
		typedef typename FEMTree< Dim , Real >::template InterpolationInfo< Real , 0 > InterpolationInfo;

		// The finite-element tracking tree node
		typedef RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type > FEMTreeNode;

		typedef typename FEMTreeInitializer< Dim , Real >::GeometryNodeType GeometryNodeType;

		// The type of the auxiliary information (including the normal)
		typedef typename std::conditional< HasAuxData , VectorTypeUnion< Real , Normal< Real , Dim > , AuxData > , Normal< Real , Dim > >::type NormalAndAuxData;

		// The type describing the sampling density
		typedef typename std::conditional< HasAuxData , ReconstructionInfo< Real , Dim , FEMSig , AuxData > , ReconstructionInfo< Real , Dim , FEMSig > >::type::DensityEstimator DensityEstimator;
		// <-- Types //
		///////////////

		// The solution info to be returned
		typename std::conditional< HasAuxData , ReconstructionInfo< Real , Dim , FEMSig , AuxData > , ReconstructionInfo< Real , Dim , FEMSig > >::type *sInfo;
		if constexpr( HasAuxData ) sInfo = new ReconstructionInfo< Real , Dim , FEMSig , AuxData >( pointStream.zero() );
		else sInfo = new ReconstructionInfo< Real , Dim , FEMSig >();

		NormalAndAuxData zeroNormalAndAuxData;
		if constexpr( HasAuxData ) zeroNormalAndAuxData = NormalAndAuxData( Normal< Real , Dim >() , sInfo->zeroAuxData );

		XForm< Real , Dim+1 > modelToUnitCube = XForm< Real , Dim+1 >::Identity();

		Profiler profiler(20);

		size_t pointCount;

		ProjectiveData< Point< Real , 2 > , Real > pointDepthAndWeight;
#ifdef SOFT_DIRICHLET
		std::vector< typename FEMTree< Dim , Real >::PointSample > *dirichletSamples = NULL;
#endif // SOFT_DIRICHLET
		DenseNodeData< GeometryNodeType , IsotropicUIntPack< Dim , FEMTrivialSignature > > geometryNodeDesignators;
		SparseNodeData< Point< Real , Dim > , NormalSigs > *normalInfo = NULL;
		std::vector< typename FEMTree< Dim , Real >::PointSample > *samples = new std::vector< typename FEMTree< Dim , Real >::PointSample >();
		std::vector< NormalAndAuxData > *sampleNormalAndAuxData = NULL;

		Real targetValue = (Real)0.5;

		// Read in the samples (and auxiliary data)
		{
			profiler.reset();

			pointStream.reset();
			sampleNormalAndAuxData = new std::vector< NormalAndAuxData >();

			modelToUnitCube = params.scale>0 ? GetPointXForm< Real , Dim >( pointStream , zeroNormalAndAuxData , params.scale ) * modelToUnitCube : modelToUnitCube;

			if( params.width>0 )
			{
				Real maxScale = 0;
				for( unsigned int i=0 ; i<Dim ; i++ ) maxScale = std::max< Real >( maxScale , (Real)1./modelToUnitCube(i,i) );
				params.depth = (unsigned int)ceil( std::max< double >( 0. , log( maxScale/params.width )/log(2.) ) );
			}
			if( params.solveDepth>params.depth )
			{
				if( params.solveDepth!=-1 ) WARN( "Solution depth cannot exceed system depth: " , params.solveDepth , " <= " , params.depth );
				params.solveDepth = params.depth;
			}
			if( params.fullDepth>params.solveDepth )
			{
				if( params.fullDepth!=-1 ) WARN( "Full depth cannot exceed system depth: " , params.fullDepth , " <= " , params.solveDepth );
				params.fullDepth = params.solveDepth;
			}
			if( params.baseDepth>params.fullDepth )
			{
				if( params.baseDepth!=-1 ) WARN( "Base depth must be smaller than full depth: " , params.baseDepth , " <= " , params.fullDepth );
				params.baseDepth = params.fullDepth;
			}
			if( params.kernelDepth==-1 ) params.kernelDepth = params.depth>2 ? params.depth-2 : 0;
			if( params.kernelDepth>params.depth )
			{
				if( params.kernelDepth!=-1 ) WARN( "Kernel depth cannot exceed system depth: " , params.kernelDepth , " <= " , params.depth );
				params.kernelDepth = params.depth;
			}

			if constexpr( HasAuxData )
			{
#ifdef DE_VIRTUALIZE_OUTPUT
				TransformedInputSampleWithDataStream< Real , Dim , AuxData , InputSampleStreamType > _pointStream( modelToUnitCube , pointStream );
#else // !DE_VIRTUALIZE_OUTPUT
				TransformedInputSampleWithDataStream< Real , Dim , AuxData > _pointStream( modelToUnitCube , pointStream );
#endif // DE_VIRTUALIZE_OUTPUT
				auto ProcessDataWithConfidence = [&]( const Point< Real , Dim > &p , NormalAndAuxData &d )
					{
						Real l = (Real)Length( d.template get<0>() );
						if( !l || !std::isfinite( l ) ) return (Real)-1.;
						return (Real)pow( l , params.confidence );
					};
				auto ProcessData = []( const Point< Real , Dim > &p , NormalAndAuxData &d )
					{
						Real l = (Real)Length( d.template get<0>() );
						if( !l || !std::isfinite( l ) ) return (Real)-1.;
						d.template get<0>() /= l;
						return (Real)1.;
					};

				typename FEMTreeInitializer< Dim , Real >::StreamInitializationData sid;
				if( params.confidence>0 ) pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< NormalAndAuxData >( sid , sInfo->tree.spaceRoot() , _pointStream , zeroNormalAndAuxData , params.depth , *samples , *sampleNormalAndAuxData , true , sInfo->tree.nodeAllocators.size() ? sInfo->tree.nodeAllocators[0] : NULL , sInfo->tree.initializer() , ProcessDataWithConfidence );
				else                      pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< NormalAndAuxData >( sid , sInfo->tree.spaceRoot() , _pointStream , zeroNormalAndAuxData , params.depth , *samples , *sampleNormalAndAuxData , true , sInfo->tree.nodeAllocators.size() ? sInfo->tree.nodeAllocators[0] : NULL , sInfo->tree.initializer() , ProcessData );
			}
			else
			{
#ifdef DE_VIRTUALIZE_OUTPUT
				TransformedInputSampleStream< Real , Dim , InputSampleStreamType > _pointStream( modelToUnitCube , pointStream );
#else // !DE_VIRTUALIZE_OUTPUT
				TransformedInputSampleStream< Real , Dim > _pointStream( modelToUnitCube , pointStream );
#endif // DE_VIRTUALIZE_OUTPUT
				auto ProcessDataWithConfidence = [&]( const Point< Real , Dim > &p , NormalAndAuxData &d )
					{
						Real l = (Real)Length( d );
						if( !l || !std::isfinite( l ) ) return (Real)-1.;
						return (Real)pow( l , params.confidence );
					};
				auto ProcessData = []( const Point< Real , Dim > &p , NormalAndAuxData &d )
					{
						Real l = (Real)Length( d );
						if( !l || !std::isfinite( l ) ) return (Real)-1.;
						d /= l;
						return (Real)1.;
					};

				typename FEMTreeInitializer< Dim , Real >::StreamInitializationData sid;
				if( params.confidence>0 ) pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< NormalAndAuxData >( sid , sInfo->tree.spaceRoot() , _pointStream , zeroNormalAndAuxData , params.depth , *samples , *sampleNormalAndAuxData , true , sInfo->tree.nodeAllocators.size() ? sInfo->tree.nodeAllocators[0] : NULL , sInfo->tree.initializer() , ProcessDataWithConfidence );
				else                      pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< NormalAndAuxData >( sid , sInfo->tree.spaceRoot() , _pointStream , zeroNormalAndAuxData , params.depth , *samples , *sampleNormalAndAuxData , true , sInfo->tree.nodeAllocators.size() ? sInfo->tree.nodeAllocators[0] : NULL , sInfo->tree.initializer() , ProcessData );
			}

			sInfo->unitCubeToModel = modelToUnitCube.inverse();

			if( params.verbose )
			{
				std::cout << "Input Points / Samples: " << pointCount << " / " << samples->size() << std::endl;
				std::cout << "# Read input into tree: " << profiler << std::endl;
			}
		}
		{
#ifdef SOFT_DIRICHLET
			InterpolationInfo *dirichletInfo = NULL;
#endif // SOFT_DIRICHLET
			DenseNodeData< Real , Sigs > constraints;
			InterpolationInfo *iInfo = NULL;
			int solveDepth = params.depth;

			sInfo->tree.resetNodeIndices( 0 , std::make_tuple() );

			// Get the kernel density estimator
			{
				profiler.reset();
				sInfo->density = sInfo->tree.template setDensityEstimator< 1 , Reconstructor::WeightDegree >( *samples , params.kernelDepth , params.samplesPerNode );
				if( params.verbose ) std::cout << "#   Got kernel density: " << profiler << std::endl;
			}

			// Transform the Hermite samples into a vector field
			{
				profiler.reset();
				normalInfo = new SparseNodeData< Point< Real , Dim > , NormalSigs >();
				std::function< bool ( NormalAndAuxData , Point< Real , Dim > & ) > ConversionFunction;
				std::function< bool ( NormalAndAuxData , Point< Real , Dim > & , Real & ) > ConversionAndBiasFunction;
				if constexpr( HasAuxData )
				{
					ConversionFunction = []( NormalAndAuxData in , Point< Real , Dim > &out )
						{
							Point< Real , Dim > n = in.template get<0>();
							Real l = (Real)Length( n );
							// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
							if( !l ) return false;
							out = n / l;
							return true;
						};
					ConversionAndBiasFunction = [&]( NormalAndAuxData in , Point< Real , Dim > &out , Real &bias )
						{
							Point< Real , Dim > n = in.template get<0>();
							Real l = (Real)Length( n );
							// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
							if( !l ) return false;
							out = n / l;
							bias = (Real)( log( l ) * params.confidenceBias / log( 1<<(Dim-1) ) );
							return true;
						};
				}
				else
				{
					// In this case NormalAndAuxData = Point< Real , Dim >
					ConversionFunction = []( NormalAndAuxData in , Point< Real , Dim > &out )
						{
							Real l = (Real)Length( in );
							// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
							if( !l ) return false;
							out = in / l;
							return true;
						};
					ConversionAndBiasFunction = [&]( NormalAndAuxData in , Point< Real , Dim > &out , Real &bias )
						{
							Real l = (Real)Length( in );
							// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
							if( !l ) return false;
							out = in / l;
							bias = (Real)( log( l ) * params.confidenceBias / log( 1<<(Dim-1) ) );
							return true;
						};
				}
				if( params.confidenceBias>0 ) *normalInfo = sInfo->tree.setInterpolatedDataField( NormalSigs() , *samples , *sampleNormalAndAuxData , sInfo->density , params.baseDepth , params.depth , params.lowDepthCutOff , pointDepthAndWeight , ConversionAndBiasFunction );
				else                          *normalInfo = sInfo->tree.setInterpolatedDataField( NormalSigs() , *samples , *sampleNormalAndAuxData , sInfo->density , params.baseDepth , params.depth , params.lowDepthCutOff , pointDepthAndWeight , ConversionFunction );
				ThreadPool::Parallel_for( 0 , normalInfo->size() , [&]( unsigned int , size_t i ){ (*normalInfo)[i] *= (Real)-1.; } );
				if( params.verbose )
				{
					std::cout << "#     Got normal field: " << profiler << std::endl;
					std::cout << "Point depth / Point weight / Estimated measure: " << pointDepthAndWeight.value()[0] << " / " << pointDepthAndWeight.value()[1] << " / " << pointCount*pointDepthAndWeight.value()[1] << std::endl;
				}
			}

			// Get the geometry designators indicating if the space node are interior to, exterior to, or contain the envelope boundary
			if( envelopeMesh )
			{
				profiler.reset();
#ifdef SOFT_DIRICHLET
				if( params.dirichletWeight>0 )
				{
					std::vector< Point< Real , Dim > > vertices( envelopeMesh->vertices.size() );
					for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i] = modelToUnitCube * envelopeMesh->vertices[i];

#if 0
					// Get the coarsest interior/boundary/exterior designators
					if( exteriorConstraintType!=ExteriorConstraint::NONE || interiorConstraintType!=InteriorConstraint::NONE )
						geometryNodeDesignators = FEMTreeInitializer< Dim , Real >::template GetGeometryNodeDesignators( sInfo->tree.spaceRoot() , vertices , envelopeMesh->simplices , params.envelopeDepth , tree.nodeAllocators , tree.initializer() );
#endif
					// Get the samples of the envelope that will act as Dirichlet point constraints and turn on the scratch flags for nodes containing constraints
					dirichletSamples = new std::vector< typename FEMTree< Dim , Real >::PointSample >();
					FEMTreeInitializer< Dim , Real >::Initialize( sInfo->tree.spaceRoot() , vertices , envelopeMesh->simplices , params.envelopeDepth , *dirichletSamples , true , sInfo->tree.nodeAllocators , sInfo->tree.initializer() );
					ThreadPool::Parallel_for( 0 , dirichletSamples->size() , [&]( unsigned int , size_t i ){ for( FEMTreeNode *node=(*dirichletSamples)[i].node ; node ; node=node->parent ) node->nodeData.setScratchFlag( true ); } );
				}
				else
#endif // SOFT_DIRICHLET
				{
					// Make the octree complete up to the base depth
					FEMTreeInitializer< Dim , Real >::Initialize( sInfo->tree.spaceRoot() , params.baseDepth , []( int , int[] ){ return true; } , sInfo->tree.nodeAllocators.size() ?  sInfo->tree.nodeAllocators[0] : NULL , sInfo->tree.initializer() );

					std::vector< Point< Real , Dim > > vertices( envelopeMesh->vertices.size() );
					for( unsigned int i=0 ; i<vertices.size() ; i++ ) vertices[i] = modelToUnitCube * envelopeMesh->vertices[i];
					geometryNodeDesignators = FEMTreeInitializer< Dim , Real >::template GetGeometryNodeDesignators( &sInfo->tree.spaceRoot() , vertices , envelopeMesh->simplices , params.baseDepth , params.envelopeDepth , sInfo->tree.nodeAllocators , sInfo->tree.initializer() );

					// Make nodes in the support of the vector field @{ExactDepth} interior
					if( params.dirichletErode )
					{
						// What to do if we find a node in the support of the vector field
						auto SetScratchFlag = [&]( FEMTreeNode *node )
							{
								if( node )
								{
									while( node->depth()>(int)params.baseDepth ) node = node->parent;
									node->nodeData.setScratchFlag( true );
								}
							};

						std::function< void ( FEMTreeNode * ) > PropagateToLeaves = [&]( const FEMTreeNode *node )
							{
								geometryNodeDesignators[ node ] = GeometryNodeType::INTERIOR;
								if( node->children ) for( int c=0 ; c<(1<<Dim) ; c++ ) PropagateToLeaves( node->children+c );
							};

						// Flags indicating if a node contains a non-zero vector field coefficient
						std::vector< bool > isVectorFieldElement( sInfo->tree.nodeCount() , false );

						// Get the set of base nodes
						std::vector< FEMTreeNode * > baseNodes;
						auto nodeFunctor = [&]( FEMTreeNode *node )
							{
								if( node->depth()==params.baseDepth ) baseNodes.push_back( node );
								return node->depth()<(int)params.baseDepth;
							};
						sInfo->tree.spaceRoot().processNodes( nodeFunctor );

						std::vector< node_index_type > vectorFieldElementCounts( baseNodes.size() );
						for( int i=0 ; i<vectorFieldElementCounts.size() ; i++ ) vectorFieldElementCounts[i] = 0;

						// In parallel, iterate over the base nodes and mark the nodes containing non-zero vector field coefficients
						ThreadPool::Parallel_for( 0 , baseNodes.size() , [&]( unsigned int t , size_t  i )
							{
								auto nodeFunctor = [&]( FEMTreeNode *node )
									{
										Point< Real , Dim > *n = (*normalInfo)( node );
										if( n && Point< Real , Dim >::SquareNorm( *n ) ) isVectorFieldElement[ node->nodeData.nodeIndex ] = true , vectorFieldElementCounts[i]++;
									};
								baseNodes[i]->processNodes( nodeFunctor );
							} );
						size_t vectorFieldElementCount = 0;
						for( int i=0 ; i<vectorFieldElementCounts.size() ; i++ ) vectorFieldElementCount += vectorFieldElementCounts[i];

						// Get the subset of nodes containing non-zero vector field coefficients and disable the "scratch" flag
						std::vector< FEMTreeNode * > vectorFieldElements;
						vectorFieldElements.reserve( vectorFieldElementCount );
						{
							std::vector< std::vector< FEMTreeNode * > > _vectorFieldElements( baseNodes.size() );
							for( int i=0 ; i<_vectorFieldElements.size() ; i++ ) _vectorFieldElements[i].reserve( vectorFieldElementCounts[i] );
							ThreadPool::Parallel_for( 0 , baseNodes.size() , [&]( unsigned int t , size_t  i )
								{
									auto nodeFunctor = [&]( FEMTreeNode *node )
										{
											if( isVectorFieldElement[ node->nodeData.nodeIndex ] ) _vectorFieldElements[i].push_back( node );
											node->nodeData.setScratchFlag( false );
										};
									baseNodes[i]->processNodes( nodeFunctor );
								} );
							for( int i=0 ; i<_vectorFieldElements.size() ; i++ ) vectorFieldElements.insert( vectorFieldElements.end() , _vectorFieldElements[i].begin() , _vectorFieldElements[i].end() );
						}

						// Set the scratch flag for the base nodes on which the vector field is supported
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] In principal, we should unlock finite elements whose support overlaps the vector field" )
#endif // SHOW_WARNINGS
						sInfo->tree.template processNeighboringLeaves< -BSplineSupportSizes< Poisson::NormalDegree >::SupportStart , BSplineSupportSizes< Poisson::NormalDegree >::SupportEnd >( &vectorFieldElements[0] , vectorFieldElements.size() , SetScratchFlag , false );

						// Set sub-trees rooted at interior nodes @ ExactDepth to interior
						ThreadPool::Parallel_for( 0 , baseNodes.size() , [&]( unsigned int , size_t  i ){ if( baseNodes[i]->nodeData.getScratchFlag() ) PropagateToLeaves( baseNodes[i] ); } );

						// Adjust the coarser node designators in case exterior nodes have become boundary.
						ThreadPool::Parallel_for( 0 , baseNodes.size() , [&]( unsigned int , size_t  i ){ FEMTreeInitializer< Dim , Real >::PullGeometryNodeDesignatorsFromFiner( baseNodes[i] , geometryNodeDesignators ); } );
						FEMTreeInitializer< Dim , Real >::PullGeometryNodeDesignatorsFromFiner( &sInfo->tree.spaceRoot() , geometryNodeDesignators , params.baseDepth );
					}
				}
				if( params.verbose ) std::cout << "#               Initialized envelope constraints: " << profiler << std::endl;
			}

			if( !params.outputDensity ){ delete sInfo->density ; sInfo->density = NULL; }
			if constexpr( HasAuxData ) sInfo->auxData = new SparseNodeData< ProjectiveData< AuxData , Real > , IsotropicUIntPack< Dim , DataSig > >( sInfo->tree.template setExtrapolatedDataField< DataSig , false , Reconstructor::WeightDegree , AuxData >( samples->size() , [&]( size_t i ) -> const typename FEMTree< Dim , Real >::PointSample & { return (*samples)[i]; } , [&]( size_t i ) -> const AuxData & { return (*sampleNormalAndAuxData)[i].template get<1>(); } , (DensityEstimator*)NULL ) );
			delete sampleNormalAndAuxData;

			// Add the interpolation constraints
			if( params.pointWeight>0 )
			{
				profiler.reset();
				if( params.exactInterpolation ) iInfo = FEMTree< Dim , Real >::template       InitializeExactPointInterpolationInfo< Real , 0 > ( sInfo->tree , *samples , Poisson::ConstraintDual< Dim , Real >( targetValue , params.pointWeight * pointDepthAndWeight.value()[1] ) , Poisson::SystemDual< Dim , Real >( params.pointWeight * pointDepthAndWeight.value()[1] ) , true , false );
				else                            iInfo = FEMTree< Dim , Real >::template InitializeApproximatePointInterpolationInfo< Real , 0 > ( sInfo->tree , *samples , Poisson::ConstraintDual< Dim , Real >( targetValue , params.pointWeight * pointDepthAndWeight.value()[1] ) , Poisson::SystemDual< Dim , Real >( params.pointWeight * pointDepthAndWeight.value()[1] ) , true , params.depth , 1 );
				if( params.verbose ) std::cout <<  "#Initialized point interpolation constraints: " << profiler << std::endl;
			}

			// Trim the tree and prepare for multigrid
			{
				profiler.reset();
				constexpr int MaxDegree = Poisson::NormalDegree > Degrees::Max() ? Poisson::NormalDegree : Degrees::Max();
				typename FEMTree< Dim , Real >::template HasNormalDataFunctor< NormalSigs > hasNormalDataFunctor( *normalInfo );
				auto hasDataFunctor = [&]( const FEMTreeNode *node ){ return hasNormalDataFunctor( node ); };
				auto addNodeFunctor = [&]( int d , const int off[Dim] ){ return d<=(int)params.fullDepth; };
				if constexpr( HasAuxData )
				{
					if( geometryNodeDesignators.size() ) sInfo->tree.template finalizeForMultigridWithDirichlet< MaxDegree , Degrees::Max() >( params.baseDepth , addNodeFunctor , hasDataFunctor , [&]( const FEMTreeNode *node ){ return node->nodeData.nodeIndex<(node_index_type)geometryNodeDesignators.size() && geometryNodeDesignators[node]==GeometryNodeType::EXTERIOR; } , std::make_tuple( iInfo ) , std::make_tuple( normalInfo , sInfo->density , sInfo->auxData , &geometryNodeDesignators ) );
					else                                 sInfo->tree.template finalizeForMultigrid             < MaxDegree , Degrees::Max() >( params.baseDepth , addNodeFunctor , hasDataFunctor ,                                                                                                                                                                                   std::make_tuple( iInfo ) , std::make_tuple( normalInfo , sInfo->density , sInfo->auxData ) );
				}
				else
				{
					if( geometryNodeDesignators.size() ) sInfo->tree.template finalizeForMultigridWithDirichlet< MaxDegree , Degrees::Max() >( params.baseDepth , addNodeFunctor , hasDataFunctor , [&]( const FEMTreeNode *node ){ return node->nodeData.nodeIndex<(node_index_type)geometryNodeDesignators.size() && geometryNodeDesignators[node]==GeometryNodeType::EXTERIOR; } , std::make_tuple( iInfo ) , std::make_tuple( normalInfo , sInfo->density , &geometryNodeDesignators ) );
					else                                 sInfo->tree.template finalizeForMultigrid             < MaxDegree , Degrees::Max() >( params.baseDepth , addNodeFunctor , hasDataFunctor ,                                                                                                                                                                                   std::make_tuple( iInfo ) , std::make_tuple( normalInfo , sInfo->density ) );
				}

				if( params.verbose ) std::cout << "#       Finalized tree: " << profiler << std::endl;
			}

			// Add the FEM constraints
			{
				profiler.reset();
				constraints = sInfo->tree.initDenseNodeData( Sigs() );

				// Add Poisson constraints
				{
					typename FEMIntegrator::template Constraint< Sigs , IsotropicUIntPack< Dim , 1 > , NormalSigs , IsotropicUIntPack< Dim , 0 > , Dim > F;
					unsigned int derivatives2[Dim];
					for( int d=0 ; d<Dim ; d++ ) derivatives2[d] = 0;
					typedef IsotropicUIntPack< Dim , 1 > Derivatives1;
					typedef IsotropicUIntPack< Dim , 0 > Derivatives2;
					for( int d=0 ; d<Dim ; d++ )
					{
						unsigned int derivatives1[Dim];
						for( int dd=0 ; dd<Dim ; dd++ ) derivatives1[dd] = dd==d ? 1 : 0;
						F.weights[d][ TensorDerivatives< Derivatives1 >::Index( derivatives1 ) ][ TensorDerivatives< Derivatives2 >::Index( derivatives2 ) ] = 1;
					}
					sInfo->tree.addFEMConstraints( F , *normalInfo , constraints , solveDepth );
				}
				if( params.verbose ) std::cout << "#  Set FEM constraints: " << profiler << std::endl;
			}

			// Free up the normal info
			delete normalInfo , normalInfo = NULL;

			if( params.pointWeight>0 )
			{
				profiler.reset();
				sInfo->tree.addInterpolationConstraints( constraints , solveDepth , std::make_tuple( iInfo ) );
				if( params.verbose ) std::cout << "#Set point constraints: " << profiler << std::endl;
			}

#ifdef SOFT_DIRICHLET
			if( dirichletSamples )
			{
				if( params.exactInterpolation ) dirichletInfo = FEMTree< Dim , Real >::template       InitializeExactPointInterpolationInfo< Real , 0 > ( sInfo->tree , *dirichletSamples , Poisson::ConstraintDual< Dim , Real >( 0. , params.dirichletWeight ) , Poisson::SystemDual< Dim , Real >( params.dirichletWeight ) , true , false );
				else                            dirichletInfo = FEMTree< Dim , Real >::template InitializeApproximatePointInterpolationInfo< Real , 0 > ( sInfo->tree , *dirichletSamples , Poisson::ConstraintDual< Dim , Real >( 0. , params.dirichletWeight ) , Poisson::SystemDual< Dim , Real >( params.dirichletWeight ) , true , 1 );
				//			if( exteriorConstraintValue ) sInfo->tree.addInterpolationConstraints( constraints , solveDepth , *dirichletInfo );
				if( params.verbose ) std::cout << "# Set exterior envelope constraints: " << profiler << std::endl;
				delete dirichletSamples , dirichletSamples = NULL;
			}
#endif // SOFT_DIRICHLET


#ifdef SOFT_DIRICHLET
			if( params.verbose )
				if( params.dirichletWeight>0 ) std::cout << "All Nodes / Active Nodes / Ghost Nodes: " << sInfo->tree.allNodes() << " / " << sInfo->tree.activeNodes() << " / " << sInfo->tree.ghostNodes() << std::endl;
				else                           std::cout << "All Nodes / Active Nodes / Ghost Nodes / Dirichlet Supported Nodes: " << sInfo->tree.allNodes() << " / " << sInfo->tree.activeNodes() << " / " << sInfo->tree.ghostNodes() << " / " << sInfo->tree.dirichletElements() << std::endl;
#else // !SOFT_DIRICHLET
			if( params.verbose ) std::cout << "All Nodes / Active Nodes / Ghost Nodes / Dirichlet Supported Nodes: " << sInfo->tree.allNodes() << " / " << sInfo->tree.activeNodes() << " / " << sInfo->tree.ghostNodes() << " / " << sInfo->tree.dirichletElements() << std::endl;
#endif // SOFT_DIRICHLET
			if( params.verbose ) std::cout << "Memory Usage: " << float( MemoryInfo::Usage())/(1<<20) << " MB" << std::endl;

			// Solve the linear system
			{
				profiler.reset();
				typename FEMTree< Dim , Real >::SolverInfo _sInfo;
				_sInfo.cgDepth = 0 , _sInfo.cascadic = true , _sInfo.vCycles = 1 , _sInfo.iters = params.iters , _sInfo.cgAccuracy = params.cgSolverAccuracy , _sInfo.verbose = params.verbose , _sInfo.showResidual = params.showResidual , _sInfo.showGlobalResidual = SHOW_GLOBAL_RESIDUAL_NONE , _sInfo.sliceBlockSize = 1;
				_sInfo.baseVCycles = params.baseVCycles;
				typename FEMIntegrator::template System< Sigs , IsotropicUIntPack< Dim , 1 > > F( { 0. , 1. } );
#ifdef SOFT_DIRICHLET
				if( dirichletInfo ) sInfo->solution = sInfo->tree.solveSystem( Sigs() , F , constraints , params.baseDepth , params.solveDepth , _sInfo , std::make_tuple( iInfo , dirichletInfo ) );
				else                sInfo->solution = sInfo->tree.solveSystem( Sigs() , F , constraints , params.baseDepth , params.solveDepth , _sInfo , std::make_tuple( iInfo ) );
#else // !SOFT_DIRICHLET
				sInfo->solution = sInfo->tree.solveSystem( Sigs() , F , constraints , params.baseDepth , params.solveDepth , _sInfo , std::make_tuple( iInfo ) );
#endif // SOFT_DIRICHLET
				if( params.verbose ) std::cout << "# Linear system solved: " << profiler << std::endl;
				if( iInfo ) delete iInfo , iInfo = NULL;
#ifdef SOFT_DIRICHLET
				if( dirichletInfo ) delete dirichletInfo , dirichletInfo = NULL;
#endif // SOFT_DIRICHLET
			}
		}

		// Get the iso-value
		{
			profiler.reset();
			double valueSum = 0 , weightSum = 0;
			typename FEMTree< Dim , Real >::template MultiThreadedEvaluator< Sigs , 0 > evaluator( &sInfo->tree , sInfo->solution );
			std::vector< double > valueSums( ThreadPool::NumThreads() , 0 ) , weightSums( ThreadPool::NumThreads() , 0 );
			ThreadPool::Parallel_for( 0 , samples->size() , [&]( unsigned int thread , size_t j )
				{
					ProjectiveData< Point< Real , Dim > , Real >& sample = (*samples)[j].sample;
					Real w = sample.weight;
					if( w>0 ) weightSums[thread] += w , valueSums[thread] += evaluator.values( sample.data / sample.weight , thread , (*samples)[j].node )[0] * w;
				} );
			for( size_t t=0 ; t<valueSums.size() ; t++ ) valueSum += valueSums[t] , weightSum += weightSums[t];
			sInfo->isoValue = (Real)( valueSum / weightSum );
			if( params.verbose )
			{
				std::cout << "Got average: " << profiler << std::endl;
				std::cout << "Iso-Value: " << sInfo->isoValue << " = " << valueSum << " / " << weightSum << std::endl;
			}
		}
		delete samples;
		return sInfo;
	}

	template< bool HasAuxData , typename Real , unsigned int Dim , unsigned int FEMSig , typename AuxData , typename InputSampleStreamType , unsigned int ... FEMSigs >
	typename std::conditional< HasAuxData , ReconstructionInfo< Real , Dim , FEMSig , AuxData > , ReconstructionInfo< Real , Dim , FEMSig > >::type *SSD::_Solve( UIntPack< FEMSigs ... >  , InputSampleStreamType &pointStream , SolutionParameters< Real > params )
	{
		static_assert( std::is_same< IsotropicUIntPack< Dim , FEMSig > , UIntPack< FEMSigs... > >::value , "[ERROR] Signatures don't match" );

		// The signature for the finite-elements representing the auxiliary data (if it's there)
		static const unsigned int DataSig = FEMDegreeAndBType< Reconstructor::DataDegree , BOUNDARY_FREE >::Signature;

		///////////////
		// Types --> //
		// The packed finite elements signature
		typedef UIntPack< FEMSigs ... > Sigs;

		// The degrees of the finite elements across the different axes
		typedef UIntPack< FEMSignature< FEMSigs >::Degree ... > Degrees;

		// The signature describing the normals elements
		typedef UIntPack< FEMDegreeAndBType< SSD::NormalDegree , DerivativeBoundary< FEMSignature< FEMSigs >::BType , 1 >::BType >::Signature ... > NormalSigs;

		// Type for tracking sample interpolation
		typedef typename FEMTree< Dim , Real >::template InterpolationInfo< Real , 1 > InterpolationInfo;

		// The finite-element tracking tree node
		typedef RegularTreeNode< Dim , FEMTreeNodeData , depth_and_offset_type > FEMTreeNode;

		// The type of the auxiliary information (including the normal)
		typedef typename std::conditional< HasAuxData , VectorTypeUnion< Real , Normal< Real , Dim > , AuxData > , Normal< Real , Dim > >::type NormalAndAuxData;

		// The type describing the sampling density
		typedef typename std::conditional< HasAuxData , ReconstructionInfo< Real , Dim , FEMSig , AuxData > , ReconstructionInfo< Real , Dim , FEMSig > >::type::DensityEstimator DensityEstimator;
		// <-- Types //
		///////////////

		// The solution info to be returned
		typename std::conditional< HasAuxData , ReconstructionInfo< Real , Dim , FEMSig , AuxData > , ReconstructionInfo< Real , Dim , FEMSig > >::type *sInfo;
		if constexpr( HasAuxData ) sInfo = new ReconstructionInfo< Real , Dim , FEMSig , AuxData >( pointStream.zero() );
		else sInfo = new ReconstructionInfo< Real , Dim , FEMSig >();

		NormalAndAuxData zeroNormalAndAuxData;
		if constexpr( HasAuxData ) zeroNormalAndAuxData = NormalAndAuxData( Normal< Real , Dim >() , sInfo->zeroAuxData );

		XForm< Real , Dim+1 > modelToUnitCube = XForm< Real , Dim+1 >::Identity();

		Profiler profiler(20);

		size_t pointCount;

		ProjectiveData< Point< Real , 2 > , Real > pointDepthAndWeight;
		SparseNodeData< Point< Real , Dim > , NormalSigs > *normalInfo = NULL;
		std::vector< typename FEMTree< Dim , Real >::PointSample > *samples = new std::vector< typename FEMTree< Dim , Real >::PointSample >();
		std::vector< NormalAndAuxData > *sampleNormalAndAuxData = NULL;

		Real targetValue = (Real)0.0;

		// Read in the samples (and auxiliary data)
		{
			profiler.reset();

			pointStream.reset();
			sampleNormalAndAuxData = new std::vector< NormalAndAuxData >();

			modelToUnitCube = params.scale>0 ? GetPointXForm< Real , Dim >( pointStream , zeroNormalAndAuxData , params.scale ) * modelToUnitCube : modelToUnitCube;

			if( params.width>0 )
			{
				Real maxScale = 0;
				for( unsigned int i=0 ; i<Dim ; i++ ) maxScale = std::max< Real >( maxScale , (Real)1./modelToUnitCube(i,i) );
				params.depth = (unsigned int)ceil( std::max< double >( 0. , log( maxScale/params.width )/log(2.) ) );
			}
			if( params.solveDepth>params.depth )
			{
				if( params.solveDepth!=-1 ) WARN( "Solution depth cannot exceed system depth: " , params.solveDepth , " <= " , params.depth );
				params.solveDepth = params.depth;
			}
			if( params.fullDepth>params.solveDepth )
			{
				if( params.fullDepth!=-1 ) WARN( "Full depth cannot exceed system depth: " , params.fullDepth , " <= " , params.solveDepth );
				params.fullDepth = params.solveDepth;
			}
			if( params.baseDepth>params.fullDepth )
			{
				if( params.baseDepth!=-1 ) WARN( "Base depth must be smaller than full depth: " , params.baseDepth , " <= " , params.fullDepth );
				params.baseDepth = params.fullDepth;
			}
			if( params.kernelDepth==-1 ) params.kernelDepth = params.depth>2 ? params.depth-2 : 0;
			if( params.kernelDepth>params.depth )
			{
				if( params.kernelDepth!=-1 ) WARN( "Kernel depth cannot exceed system depth: " , params.kernelDepth , " <= " , params.depth );
				params.kernelDepth = params.depth;
			}

			if constexpr( HasAuxData )
			{
#ifdef DE_VIRTUALIZE_OUTPUT
				TransformedInputSampleWithDataStream< Real , Dim , AuxData , InputSampleStreamType > _pointStream( modelToUnitCube , pointStream );
#else // !DE_VIRTUALIZE_OUTPUT
				TransformedInputSampleWithDataStream< Real , Dim , AuxData > _pointStream( modelToUnitCube , pointStream );
#endif // DE_VIRTUALIZE_OUTPUT
				auto ProcessDataWithConfidence = [&]( const Point< Real , Dim > &p , NormalAndAuxData &d )
					{
						Real l = (Real)Length( d.template get<0>() );
						if( !l || !std::isfinite( l ) ) return (Real)-1.;
						return (Real)pow( l , params.confidence );
					};
				auto ProcessData = []( const Point< Real , Dim > &p , NormalAndAuxData &d )
					{
						Real l = (Real)Length( d.template get<0>() );
						if( !l || !std::isfinite( l ) ) return (Real)-1.;
						d.template get<0>() /= l;
						return (Real)1.;
					};

				typename FEMTreeInitializer< Dim , Real >::StreamInitializationData sid;
				if( params.confidence>0 ) pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< NormalAndAuxData >( sid , sInfo->tree.spaceRoot() , _pointStream , zeroNormalAndAuxData , params.depth , *samples , *sampleNormalAndAuxData , true , sInfo->tree.nodeAllocators.size() ? sInfo->tree.nodeAllocators[0] : NULL , sInfo->tree.initializer() , ProcessDataWithConfidence );
				else                      pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< NormalAndAuxData >( sid , sInfo->tree.spaceRoot() , _pointStream , zeroNormalAndAuxData , params.depth , *samples , *sampleNormalAndAuxData , true , sInfo->tree.nodeAllocators.size() ? sInfo->tree.nodeAllocators[0] : NULL , sInfo->tree.initializer() , ProcessData );
			}
			else
			{
#ifdef DE_VIRTUALIZE_OUTPUT
				TransformedInputSampleStream< Real , Dim , InputSampleStreamType > _pointStream( modelToUnitCube , pointStream );
#else // !DE_VIRTUALIZE_OUTPUT
				TransformedInputSampleStream< Real , Dim > _pointStream( modelToUnitCube , pointStream );
#endif // DE_VIRTUALIZE_OUTPUT
				auto ProcessDataWithConfidence = [&]( const Point< Real , Dim > &p , NormalAndAuxData &d )
					{
						Real l = (Real)Length( d );
						if( !l || !std::isfinite( l ) ) return (Real)-1.;
						return (Real)pow( l , params.confidence );
					};
				auto ProcessData = []( const Point< Real , Dim > &p , NormalAndAuxData &d )
					{
						Real l = (Real)Length( d );
						if( !l || !std::isfinite( l ) ) return (Real)-1.;
						d /= l;
						return (Real)1.;
					};

				typename FEMTreeInitializer< Dim , Real >::StreamInitializationData sid;
				if( params.confidence>0 ) pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< NormalAndAuxData >( sid , sInfo->tree.spaceRoot() , _pointStream , zeroNormalAndAuxData , params.depth , *samples , *sampleNormalAndAuxData , true , sInfo->tree.nodeAllocators.size() ? sInfo->tree.nodeAllocators[0] : NULL , sInfo->tree.initializer() , ProcessDataWithConfidence );
				else                      pointCount = FEMTreeInitializer< Dim , Real >::template Initialize< NormalAndAuxData >( sid , sInfo->tree.spaceRoot() , _pointStream , zeroNormalAndAuxData , params.depth , *samples , *sampleNormalAndAuxData , true , sInfo->tree.nodeAllocators.size() ? sInfo->tree.nodeAllocators[0] : NULL , sInfo->tree.initializer() , ProcessData );
			}

			sInfo->unitCubeToModel = modelToUnitCube.inverse();

			if( params.verbose )
			{
				std::cout << "Input Points / Samples: " << pointCount << " / " << samples->size() << std::endl;
				std::cout << "# Read input into tree: " << profiler << std::endl;
			}
		}
		{
			DenseNodeData< Real , Sigs > constraints;
			InterpolationInfo *iInfo = NULL;
			int solveDepth = params.depth;

			sInfo->tree.resetNodeIndices( 0 , std::make_tuple() );

			// Get the kernel density estimator
			{
				profiler.reset();
				sInfo->density = sInfo->tree.template setDensityEstimator< 1 , Reconstructor::WeightDegree >( *samples , params.kernelDepth , params.samplesPerNode );
				if( params.verbose ) std::cout << "#   Got kernel density: " << profiler << std::endl;
			}

			// Transform the Hermite samples into a vector field
			{
				profiler.reset();
				normalInfo = new SparseNodeData< Point< Real , Dim > , NormalSigs >();
				std::function< bool ( NormalAndAuxData , Point< Real , Dim > & ) > ConversionFunction;
				std::function< bool ( NormalAndAuxData , Point< Real , Dim > & , Real & ) > ConversionAndBiasFunction;
				if constexpr( HasAuxData )
				{
					ConversionFunction = []( NormalAndAuxData in , Point< Real , Dim > &out )
						{
							Point< Real , Dim > n = in.template get<0>();
							Real l = (Real)Length( n );
							// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
							if( !l ) return false;
							out = n / l;
							return true;
						};
					ConversionAndBiasFunction = [&]( NormalAndAuxData in , Point< Real , Dim > &out , Real &bias )
						{
							Point< Real , Dim > n = in.template get<0>();
							Real l = (Real)Length( n );
							// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
							if( !l ) return false;
							out = n / l;
							bias = (Real)( log( l ) * params.confidenceBias / log( 1<<(Dim-1) ) );
							return true;
						};
				}
				else
				{
					// In this case NormalAndAuxData = Point< Real , Dim >
					ConversionFunction = []( NormalAndAuxData in , Point< Real , Dim > &out )
						{
							Real l = (Real)Length( in );
							// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
							if( !l ) return false;
							out = in / l;
							return true;
						};
					ConversionAndBiasFunction = [&]( NormalAndAuxData in , Point< Real , Dim > &out , Real &bias )
						{
							Real l = (Real)Length( in );
							// It is possible that the samples have non-zero normals but there are two co-located samples with negative normals...
							if( !l ) return false;
							out = in / l;
							bias = (Real)( log( l ) * params.confidenceBias / log( 1<<(Dim-1) ) );
							return true;
						};
				}
				if( params.confidenceBias>0 ) *normalInfo = sInfo->tree.setInterpolatedDataField( NormalSigs() , *samples , *sampleNormalAndAuxData , sInfo->density , params.baseDepth , params.depth , params.lowDepthCutOff , pointDepthAndWeight , ConversionAndBiasFunction );
				else                          *normalInfo = sInfo->tree.setInterpolatedDataField( NormalSigs() , *samples , *sampleNormalAndAuxData , sInfo->density , params.baseDepth , params.depth , params.lowDepthCutOff , pointDepthAndWeight , ConversionFunction );
				if( params.verbose )
				{
					std::cout << "#     Got normal field: " << profiler << std::endl;
					std::cout << "Point depth / Point weight / Estimated measure: " << pointDepthAndWeight.value()[0] << " / " << pointDepthAndWeight.value()[1] << " / " << pointCount*pointDepthAndWeight.value()[1] << std::endl;
				}
			}

			if( !params.outputDensity ){ delete sInfo->density ; sInfo->density = NULL; }
			if constexpr( HasAuxData ) sInfo->auxData = new SparseNodeData< ProjectiveData< AuxData , Real > , IsotropicUIntPack< Dim , DataSig > >( sInfo->tree.template setExtrapolatedDataField< DataSig , false , Reconstructor::WeightDegree , AuxData >( samples->size() , [&]( size_t i ) -> const typename FEMTree< Dim , Real >::PointSample & { return (*samples)[i]; } , [&]( size_t i ) -> const AuxData & { return (*sampleNormalAndAuxData)[i].template get<1>(); } , (DensityEstimator*)NULL ) );

				// Add the interpolation constraints
				if( params.pointWeight>0 || params.gradientWeight>0 )
				{
					profiler.reset();
						if constexpr( HasAuxData )
						{
							if( params.exactInterpolation ) iInfo = FEMTree< Dim , Real >::template       InitializeExactPointAndDataInterpolationInfo< Real , NormalAndAuxData , 1 >( sInfo->tree , *samples , GetPointer( *sampleNormalAndAuxData ) , SSD::ConstraintDual< Dim , Real , AuxData >( targetValue , params.pointWeight * pointDepthAndWeight.value()[1] , params.gradientWeight * pointDepthAndWeight.value()[1]  ) , SSD::SystemDual< Dim , Real , AuxData >( params.pointWeight * pointDepthAndWeight.value()[1] , params.gradientWeight * pointDepthAndWeight.value()[1] ) , true , false );
							else                            iInfo = FEMTree< Dim , Real >::template InitializeApproximatePointAndDataInterpolationInfo< Real , NormalAndAuxData , 1 >( sInfo->tree , *samples , GetPointer( *sampleNormalAndAuxData ) , SSD::ConstraintDual< Dim , Real , AuxData >( targetValue , params.pointWeight * pointDepthAndWeight.value()[1] , params.gradientWeight * pointDepthAndWeight.value()[1]  ) , SSD::SystemDual< Dim , Real , AuxData >( params.pointWeight * pointDepthAndWeight.value()[1] , params.gradientWeight * pointDepthAndWeight.value()[1] ) , true , params.depth , 1 );
						}
						else
						{
							if( params.exactInterpolation ) iInfo = FEMTree< Dim , Real >::template       InitializeExactPointAndDataInterpolationInfo< Real , Point< Real , Dim > , 1 >( sInfo->tree , *samples , GetPointer( *sampleNormalAndAuxData ) , SSD::ConstraintDual< Dim , Real >( targetValue , params.pointWeight * pointDepthAndWeight.value()[1] , params.gradientWeight * pointDepthAndWeight.value()[1]  ) , SSD::SystemDual< Dim , Real >( params.pointWeight * pointDepthAndWeight.value()[1] , params.gradientWeight * pointDepthAndWeight.value()[1] ) , true , false );
							else                            iInfo = FEMTree< Dim , Real >::template InitializeApproximatePointAndDataInterpolationInfo< Real , Point< Real , Dim > , 1 >( sInfo->tree , *samples , GetPointer( *sampleNormalAndAuxData ) , SSD::ConstraintDual< Dim , Real >( targetValue , params.pointWeight * pointDepthAndWeight.value()[1] , params.gradientWeight * pointDepthAndWeight.value()[1]  ) , SSD::SystemDual< Dim , Real >( params.pointWeight * pointDepthAndWeight.value()[1] , params.gradientWeight * pointDepthAndWeight.value()[1] ) , true , params.depth , 1 );
						}
					if( params.verbose ) std::cout <<  "#Initialized point interpolation constraints: " << profiler << std::endl;
				}

			delete sampleNormalAndAuxData;


			// Trim the tree and prepare for multigrid
			{
				profiler.reset();
				constexpr int MaxDegree = SSD::NormalDegree > Degrees::Max() ? SSD::NormalDegree : Degrees::Max();
				typename FEMTree< Dim , Real >::template HasNormalDataFunctor< NormalSigs > hasNormalDataFunctor( *normalInfo );
				auto hasDataFunctor = [&]( const FEMTreeNode *node ){ return hasNormalDataFunctor( node ); };
				auto addNodeFunctor = [&]( int d , const int off[Dim] ){ return d<=(int)params.fullDepth; };
				if constexpr( HasAuxData ) sInfo->tree.template finalizeForMultigrid< MaxDegree , Degrees::Max() >( params.baseDepth , addNodeFunctor , hasDataFunctor , std::make_tuple( iInfo ) , std::make_tuple( normalInfo , sInfo->density , sInfo->auxData ) );
				else sInfo->tree.template finalizeForMultigrid< MaxDegree , Degrees::Max() >( params.baseDepth , addNodeFunctor , hasDataFunctor , std::make_tuple( iInfo ) , std::make_tuple( normalInfo , sInfo->density ) );

				if( params.verbose ) std::cout << "#       Finalized tree: " << profiler << std::endl;
			}

			// Free up the normal info
			delete normalInfo , normalInfo = NULL;

			if( params.pointWeight>0 || params.gradientWeight>0 )
			{
				profiler.reset();
				constraints = sInfo->tree.initDenseNodeData( Sigs() );
				sInfo->tree.addInterpolationConstraints( constraints , solveDepth , std::make_tuple( iInfo ) );
				if( params.verbose ) std::cout << "#Set point constraints: " << profiler << std::endl;
			}

			if( params.verbose ) std::cout << "All Nodes / Active Nodes / Ghost Nodes: " << sInfo->tree.allNodes() << " / " << sInfo->tree.activeNodes() << " / " << sInfo->tree.ghostNodes() << std::endl;
			if( params.verbose ) std::cout << "Memory Usage: " << float( MemoryInfo::Usage())/(1<<20) << " MB" << std::endl;

			// Solve the linear system
			{
				profiler.reset();
				typename FEMTree< Dim , Real >::SolverInfo _sInfo;
				_sInfo.cgDepth = 0 , _sInfo.cascadic = true , _sInfo.vCycles = 1 , _sInfo.iters = params.iters , _sInfo.cgAccuracy = params.cgSolverAccuracy , _sInfo.verbose = params.verbose , _sInfo.showResidual = params.showResidual , _sInfo.showGlobalResidual = SHOW_GLOBAL_RESIDUAL_NONE , _sInfo.sliceBlockSize = 1;
				_sInfo.baseVCycles = params.baseVCycles;
				typename FEMIntegrator::template System< Sigs , IsotropicUIntPack< Dim , 2 > > F( { 0. , 0. , (double)params.biLapWeight } );
				sInfo->solution = sInfo->tree.solveSystem( Sigs() , F , constraints , params.baseDepth , params.solveDepth , _sInfo , std::make_tuple( iInfo ) );
				if( params.verbose ) std::cout << "# Linear system solved: " << profiler << std::endl;
				if( iInfo ) delete iInfo , iInfo = NULL;
			}
		}

		// Get the iso-value
		{
			profiler.reset();
			double valueSum = 0 , weightSum = 0;
			typename FEMTree< Dim , Real >::template MultiThreadedEvaluator< Sigs , 0 > evaluator( &sInfo->tree , sInfo->solution );
			std::vector< double > valueSums( ThreadPool::NumThreads() , 0 ) , weightSums( ThreadPool::NumThreads() , 0 );
			ThreadPool::Parallel_for( 0 , samples->size() , [&]( unsigned int thread , size_t j )
				{
					ProjectiveData< Point< Real , Dim > , Real >& sample = (*samples)[j].sample;
					Real w = sample.weight;
					if( w>0 ) weightSums[thread] += w , valueSums[thread] += evaluator.values( sample.data / sample.weight , thread , (*samples)[j].node )[0] * w;
				} );
			for( size_t t=0 ; t<valueSums.size() ; t++ ) valueSum += valueSums[t] , weightSum += weightSums[t];
			sInfo->isoValue = (Real)( valueSum / weightSum );
			if( params.verbose )
			{
				std::cout << "Got average: " << profiler << std::endl;
				std::cout << "Iso-Value: " << sInfo->isoValue << " = " << valueSum << " / " << weightSum << std::endl;
			}
		}
		delete samples;
		return sInfo;
	}
}


#endif // RECONSTRUCTORS_INCLUDED