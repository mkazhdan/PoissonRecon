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

#ifndef RECONSTRUCTORS_STREAMS_INCLUDED
#define RECONSTRUCTORS_STREAMS_INCLUDED

#include <mutex>
#include "FEMTree.h"
#include "MyExceptions.h"
#include "Array.h"

namespace PoissonRecon
{
	namespace Reconstructor
	{
		template< typename Real , unsigned int Dim > using Position = Point< Real , Dim >;
		template< typename Real , unsigned int Dim > using Normal   = Point< Real , Dim >;
		template< typename Real , unsigned int Dim > using Gradient = Point< Real , Dim >;

		// Basic types
		template< typename Real , unsigned int Dim , typename ... Other >
		using Sample = std::conditional_t< sizeof...(Other)==0 , VectorTypeUnion< Real , Position< Real , Dim > , Normal< Real , Dim > > , VectorTypeUnion< Real , Position< Real , Dim > , VectorTypeUnion< Real , Normal< Real , Dim > , Other... > > >;

		// The types output by the extractor
		template< typename Real , unsigned int Dim , typename ... Other >
		using LevelSetVertex = typename LevelSetExtractor< Real , Dim , Other... >::Vertex;
		template< typename Real , unsigned int Dim , typename ... Other >
		using LevelSetIndexedVertex = typename LevelSetExtractor< Real , Dim , Other... >::IndexedVertex;

		// Stream types
		template< typename Real , unsigned int Dim , typename ... Other > struct InputSampleStream;
		template< typename Real , unsigned int Dim , typename InputStream , typename ... Other > struct TransformedInputSampleStream;

		template< typename Real , unsigned int Dim , typename ... Other > struct OutputLevelSetVertexStream;
		template< typename Real , unsigned int Dim , typename ... Other > struct OutputIndexedLevelSetVertexStream;
		template< typename Real , unsigned int Dim , typename ... Other > struct TransformedOutputLevelSetVertexStream;
		template< typename Real , unsigned int Dim , typename ... Other > struct TransformedOutputIndexedLevelSetVertexStream;
		template< typename Vertex , typename Real , unsigned int Dim , typename ... Other > struct OutputLevelSetVertexStreamWrapper;
		template< typename Vertex , typename Real , unsigned int Dim , typename ... Other > struct OutputIndexedLevelSetVertexStreamWrapper;

		template< typename Real , unsigned int Dim , bool HasGradients , bool HasDensity , typename ... Other > struct OutputLevelSetVertexInfo;
		template< typename Real , unsigned int Dim , bool HasGradients , bool HasDensity , typename ... Other > struct OutputIndexedLevelSetVertexInfo;

		template< unsigned int FaceDim , typename T=node_index_type > using Face = std::conditional_t< FaceDim==2 , std::vector< T > , std::conditional_t< FaceDim==1 , std::pair< T , T > , void * > >;
		template< unsigned int FaceDim > using  InputFaceStream =  InputDataStream< Face< FaceDim > >;
		template< unsigned int FaceDim > using OutputFaceStream = OutputDataStream< Face< FaceDim > >;

		/////////////////////////////////////
		// Value Interpolation Data Stream //
		/////////////////////////////////////
		template< typename Real , unsigned int Dim >
		struct ValueInterpolationStream : public InputDataStream< VectorTypeUnion< Real , Point< Real , Dim > , Real > >
		{
			// Functionality to reset the stream to the start
			virtual void reset( void ) = 0;

			// Functionality to extract the next position/normal pair.
			// The method returns true if there was another point in the stream to read, and false otherwise
			virtual bool base_read(                       Position< Real , Dim > &p , Real &v ) = 0;
			virtual bool base_read( unsigned int thread , Position< Real , Dim > &p , Real &v ) = 0;
			// Implementation of InputDataStream::read
			bool base_read(                       VectorTypeUnion< Real , Point< Real , Dim > , Real > &s ){ return base_read(          s.template get<0>() , s.template get<1>() ); }
			bool base_read( unsigned int thread , VectorTypeUnion< Real , Point< Real , Dim > , Real > &s ){ return base_read( thread , s.template get<0>() , s.template get<1>() ); }
		};

		/////////////////////////////////////////////////
		// Transformed Value Interpolation Data Stream //
		/////////////////////////////////////////////////
		template< typename Real , unsigned int Dim , typename InputStream >
		struct TransformedValueInterpolationStream : public ValueInterpolationStream< Real , Dim >
		{
			static_assert( std::is_base_of< ValueInterpolationStream< Real , Dim > , InputStream >::value , "[ERROR] Unexpected stream type" );

			// A constructor initialized with the transformation to be applied to the samples, and a sample stream
			TransformedValueInterpolationStream( XForm< Real , Dim+1 > xForm , InputStream &stream ) : _stream(stream) , _xForm(xForm) {}

			// Functionality to reset the stream to the start
			void reset( void ){ _stream.reset(); }

			// Functionality to extract the next position/normal pair.
			// The method returns true if there was another point in the stream to read, and false otherwise
			bool base_read( Position< Real , Dim > &p , Real &v )
			{
				VectorTypeUnion< Real , Point< Real , Dim > , Real > s;
				bool ret = _stream.read( s );
				if( ret ) p = _xForm * s.template get<0>() , v = s.template get<1>();
				return ret;
			}
			bool base_read( unsigned int thread , Position< Real , Dim > &p , Real &v )
			{
				VectorTypeUnion< Real , Point< Real , Dim > , Real > s;
				bool ret = _stream.read( thread , s );
				if( ret ) p = _xForm * s.template get<0>() , v = s.template get<1>();
				return ret;
			}

		protected:
			// A reference to the underlying stream
			InputStream &_stream;

			// The affine transformation to be applied to the positions
			XForm< Real , Dim+1 > _xForm;
		};

		///////////////////////////
		// Oriented Point Stream //
		///////////////////////////
		template< typename Real , unsigned int Dim >
		struct InputSampleStream< Real , Dim > : public InputDataStream< Sample< Real , Dim > >
		{
			// Functionality to reset the stream to the start
			virtual void reset( void ) = 0;

			// Functionality to extract the next position/normal pair.
			// The method returns true if there was another point in the stream to read, and false otherwise
			virtual bool base_read(                       Position< Real , Dim > &p , Normal< Real , Dim > &n ) = 0;
			virtual bool base_read( unsigned int thread , Position< Real , Dim > &p , Normal< Real , Dim > &n ) = 0;
			// Implementation of InputDataStream::read
			bool base_read(                       Sample< Real , Dim > &s ){ return base_read(          s.template get<0>() , s.template get<1>() ); }
			bool base_read( unsigned int thread , Sample< Real , Dim > &s ){ return base_read( thread , s.template get<0>() , s.template get<1>() ); }
		};

		///////////////////////////////////
		// Oriented Point w/ Data Stream //
		///////////////////////////////////
		template< typename Real , unsigned int Dim , typename Data >
		struct InputSampleStream< Real , Dim , Data > : public InputDataStream< Sample< Real , Dim  , Data > >			
		{
			// Functionality to reset the stream to the start
			virtual void reset( void ) = 0;

			// Functionality to extract the next position/normal pair.
			// The method returns true if there was another point in the stream to read, and false otherwise
			virtual bool base_read(                       Position< Real , Dim > &p , Normal< Real , Dim > &n , Data &d ) = 0;
			virtual bool base_read( unsigned int thread , Position< Real , Dim > &p , Normal< Real , Dim > &n , Data &d ) = 0;
			bool base_read(                       Sample< Real , Dim , Data > &s ){ return base_read(          s.template get<0>() , s.template get<1>().template get<0>() , s.template get<1>().template get<1>() ); }
			bool base_read( unsigned int thread , Sample< Real , Dim , Data > &s ){ return base_read( thread , s.template get<0>() , s.template get<1>().template get<0>() , s.template get<1>().template get<1>() ); }
		};


		///////////////////////////////////////
		// Transformed Oriented Point Stream //
		///////////////////////////////////////
		template< typename Real , unsigned int Dim , typename InputStream >
		struct TransformedInputSampleStream< Real , Dim , InputStream > : public InputSampleStream< Real , Dim >
		{
			static_assert( std::is_base_of< InputSampleStream< Real , Dim > , InputStream >::value , "[ERROR] Unexpected stream type" );
			// A constructor initialized with the transformation to be applied to the samples, and a sample stream
			TransformedInputSampleStream( XForm< Real , Dim+1 > xForm , InputStream &stream ) : _stream(stream) , _positionXForm(xForm)
			{
				_normalXForm = XForm< Real , Dim > ( xForm ).inverse().transpose() * (Real)pow( fabs( xForm.determinant() ) , 1./Dim );
			}

			// Functionality to reset the stream to the start
			void reset( void ){ _stream.reset(); }

			// Functionality to extract the next position/normal pair.
			// The method returns true if there was another point in the stream to read, and false otherwise
			bool base_read( Position< Real , Dim > &p , Normal< Real , Dim > &n )
			{
				Sample< Real , Dim > s;
				bool ret = _stream.read( s );
				if( ret ) p = _positionXForm * s.template get<0>() , n = _normalXForm * s.template get<1>();
				return ret;
			}
			bool base_read( unsigned int thread , Position< Real , Dim > &p , Normal< Real , Dim > &n )
			{
				Sample< Real , Dim > s;
				bool ret = _stream.read( thread , s );
				if( ret ) p = _positionXForm * s.template get<0>() , n = _normalXForm * s.template get<1>();
				return ret;
			}

		protected:
			// A reference to the underlying stream
			InputStream &_stream;

			// The affine transformation to be applied to the positions
			XForm< Real , Dim+1 > _positionXForm;

			// The linear transformation to be applied to the normals
			XForm< Real , Dim > _normalXForm;
		};

		///////////////////////////////////////////////
		// Transformed Oriented Point w/ Data Stream //
		///////////////////////////////////////////////
		template< typename Real , unsigned int Dim , typename InputStream , typename Data >
		struct TransformedInputSampleStream< Real , Dim , InputStream , Data > : public InputSampleStream< Real , Dim , Data >
		{
			static_assert( std::is_base_of< InputSampleStream< Real , Dim , Data > , InputStream >::value , "[ERROR] Unexpected stream type" );

			// A constructor initialized with an instance of "zero" data
			TransformedInputSampleStream( XForm< Real , Dim+1 > xForm , InputStream &stream ) : _stream(stream) , _positionXForm(xForm)
			{
				_normalXForm = XForm< Real , Dim > ( xForm ).inverse().transpose() * (Real)pow( xForm.determinant() , 1./Dim );
			}

			// Functionality to reset the stream to the start
			void reset( void ){ _stream.reset(); }

			// Functionality to extract the next position/normal pair.
			// The method returns true if there was another point in the stream to read, and false otherwise
			bool base_read( Position< Real , Dim > &p , Normal< Real , Dim > &n , Data &d )
			{
				Sample< Real , Dim , Data > s( Position< Real , Dim >() , VectorTypeUnion< Real , Normal< Real , Dim > , Data  >( Normal< Real , Dim >() , d ) );
				bool ret = _stream.read( s );
				if( ret ) p = _positionXForm * s.template get<0>() , n = _normalXForm * s.template get<1>().template get<0>() , d = s.template get<1>().template get<1>();
				return ret;
			}
			bool base_read( unsigned int thread , Position< Real , Dim > &p , Normal< Real , Dim > &n , Data &d )
			{
				Sample< Real , Dim , Data > s( Position< Real , Dim >() , VectorTypeUnion< Real , Normal< Real , Dim > , Data  >( Normal< Real , Dim >() , d ) );
				bool ret = _stream.read( thread , s );
				if( ret ) p = _positionXForm * s.template get<0>() , n = _normalXForm * s.template get<1>().template get<0>() , d = s.template get<1>().template get<1>();
				return ret;
			}

		protected:
			// A reference to the underlying stream
			InputStream &_stream;

			// The affine transformation to be applied to the positions
			XForm< Real , Dim+1 > _positionXForm;

			// The linear transformation to be applied to the normals
			XForm< Real , Dim > _normalXForm;
		};


		//////////////////////////////////////////////////////////////////////////////
		// Vertex Stream:                                                           //
		// Looks like a LevelSetVertexStream, but supports component-wise insertion //
		//////////////////////////////////////////////////////////////////////////////
		template< typename Real , unsigned int Dim >
		struct OutputLevelSetVertexStream< Real , Dim > : public OutputDataStream< LevelSetVertex< Real , Dim > >			
		{
			// Need to provide access to base write for counter support
			using OutputDataStream< LevelSetVertex< Real , Dim > >::write;

			// Functionality to insert the next vertex
			virtual void base_write(                       Position< Real , Dim > p , Gradient< Real , Dim > g , Real w ) = 0;
			virtual void base_write( unsigned int thread , Position< Real , Dim > p , Gradient< Real , Dim > g , Real w )
			{
				std::lock_guard< std::mutex > guard(_insertionMutex);
				base_write( p , g , w );
			}
			void base_write(                       const LevelSetVertex< Real , Dim > &v ){ base_write(          v.template get<0>() , v.template get<1>() , v.template get<2>() ); }
			void base_write( unsigned int thread , const LevelSetVertex< Real , Dim > &v ){ base_write( thread , v.template get<0>() , v.template get<1>() , v.template get<2>() ); }
		protected:
			std::mutex _insertionMutex;
		};

		///////////////////////////
		// Indexed Vertex Stream //
		///////////////////////////
		template< typename Real , unsigned int Dim >
		struct OutputIndexedLevelSetVertexStream< Real , Dim > : public OutputDataStream< LevelSetIndexedVertex< Real , Dim > >
		{
			// Need to provide access to base write for counter support
			using OutputDataStream< LevelSetIndexedVertex< Real , Dim > >::write;

			// Functionality to insert the next vertex
			virtual void base_write(                       node_index_type idx , Position< Real , Dim > p , Gradient< Real , Dim > g , Real w ) = 0;
			virtual void base_write( unsigned int thread , node_index_type idx , Position< Real , Dim > p , Gradient< Real , Dim > g , Real w )
			{
				std::lock_guard< std::mutex > guard(_insertionMutex);
				base_write( idx , p , g , w );
			}
			void base_write( const LevelSetIndexedVertex< Real , Dim > &v )
			{
				base_write( v.first , v.second.template get<0>() , v.second.template get<1>() , v.second.template get<2>() );
			}
			void base_write( unsigned int thread , const LevelSetIndexedVertex< Real , Dim > &v )
			{
				base_write( thread , v.first , v.second.template get<0>() , v.second.template get<1>() , v.second.template get<2>() );
			}
		protected:
			std::mutex _insertionMutex;
		};

		///////////////////////////
		// Vertex w/ Data Stream //
		///////////////////////////
		template< typename Real , unsigned int Dim , typename Data >
		struct OutputLevelSetVertexStream< Real , Dim , Data > : public OutputDataStream< LevelSetVertex< Real , Dim , Data > >
		{
			// Need to provide access to base write for counter support
			using OutputDataStream< LevelSetVertex< Real , Dim , Data > >::write;

			// Functionality to insert the next vertex
			virtual void base_write(                       Position< Real , Dim > p , Gradient< Real , Dim > g , Real w , Data d ) = 0;
			virtual void base_write( unsigned int thread , Position< Real , Dim > p , Gradient< Real , Dim > g , Real w , Data d )
			{
				std::lock_guard< std::mutex > guard( _insertionMutex );
				base_write( p , g , w , d );
			}
			void base_write(                       const LevelSetVertex< Real , Dim , Data > &v ){ return base_write(          v.template get<0>() , v.template get<1>() , v.template get<2>() , v.template get<3>() ); }
			void base_write( unsigned int thread , const LevelSetVertex< Real , Dim , Data > &v ){ return base_write( thread , v.template get<0>() , v.template get<1>() , v.template get<2>() , v.template get<3>() ); }
		protected:
			std::mutex _insertionMutex;
		};

		///////////////////////////////////
		// Indexed vertex w/ Data Stream //
		///////////////////////////////////
		template< typename Real , unsigned int Dim , typename Data >
		struct OutputIndexedLevelSetVertexStream< Real , Dim , Data > : public OutputDataStream< LevelSetIndexedVertex< Real , Dim , Data > >
		{
			// Need to provide access to base write for counter support
			using OutputDataStream< LevelSetIndexedVertex< Real , Dim , Data > >::write;

			// Functionality to insert the next vertex
			virtual void base_write(                       node_index_type idx , Position< Real , Dim > p , Gradient< Real , Dim > g , Real w , Data d ) = 0;
			virtual void base_write( unsigned int thread , node_index_type idx , Position< Real , Dim > p , Gradient< Real , Dim > g , Real w , Data d )
			{
				std::lock_guard< std::mutex > guard( _insertionMutex );
				base_write( idx , p , g , w , d );
			}
			void base_write(                       const LevelSetIndexedVertex< Real , Dim , Data > &v ){ return base_write(          v.first , v.second.template get<0>() , v.second.template get<1>() , v.second.template get<2>() , v.second.template get<3>() ); }
			void base_write( unsigned int thread , const LevelSetIndexedVertex< Real , Dim , Data > &v ){ return base_write( thread , v.first , v.second.template get<0>() , v.second.template get<1>() , v.second.template get<2>() , v.second.template get<3>() ); }
		protected:
			std::mutex _insertionMutex;
		};

		///////////////////////////////
		// Transformed Vertex Stream //
		///////////////////////////////
		template< typename Real , unsigned int Dim >
		struct TransformedOutputLevelSetVertexStream< Real , Dim > : public OutputLevelSetVertexStream< Real , Dim >
		{
			// A constructor initialized with the transformation to be applied to the samples, and a sample stream
			TransformedOutputLevelSetVertexStream( XForm< Real , Dim+1 > xForm , OutputLevelSetVertexStream< Real , Dim > &stream ) : _stream(stream) , _positionXForm(xForm)
			{
				_gradientXForm = XForm< Real , Dim > ( xForm ).inverse().transpose() * (Real)pow( xForm.determinant() , 1./Dim );
			}

			// Need to write the union to ensure that the counter gets set
			void base_write(                       Position< Real , Dim > p , Normal< Real , Dim > g , Real w ){ _stream.write(          LevelSetVertex< Real , Dim >( _positionXForm * p , _gradientXForm * g , w ) ); }
			void base_write( unsigned int thread , Position< Real , Dim > p , Normal< Real , Dim > g , Real w ){ _stream.write( thread , LevelSetVertex< Real , Dim >( _positionXForm * p , _gradientXForm * g , w ) ); }

		protected:
			// A reference to the underlying stream
			OutputLevelSetVertexStream< Real , Dim > &_stream;

			// The affine transformation to be applied to the positions
			XForm< Real , Dim+1 > _positionXForm;

			// The linear transformation to be applied to the normals
			XForm< Real , Dim > _gradientXForm;
		};


		///////////////////////////////////////
		// Transformed Vertex w/ Data Stream //
		///////////////////////////////////////
		template< typename Real , unsigned int Dim , typename Data >
		struct TransformedOutputLevelSetVertexStream< Real , Dim , Data > : public OutputLevelSetVertexStream< Real , Dim , Data >
		{
			// A constructor initialized with the transformation to be applied to the samples, and a sample stream
			TransformedOutputLevelSetVertexStream( XForm< Real , Dim+1 > xForm , OutputLevelSetVertexStream< Real , Dim , Data > &stream ) : _stream(stream) , _positionXForm(xForm)
			{
				_gradientXForm = XForm< Real , Dim > ( xForm ).inverse().transpose() * (Real)pow( xForm.determinant() , 1./Dim );
			}

			void base_write(                       Position< Real , Dim > p , Normal< Real , Dim > g , Real w , Data d ){ _stream.write(          LevelSetVertex< Real , Dim , Data >( _positionXForm * p , _gradientXForm * g , w , d ) ); }
			void base_write( unsigned int thread , Position< Real , Dim > p , Normal< Real , Dim > g , Real w , Data d ){ _stream.write( thread , LevelSetVertex< Real , Dim , Data >( _positionXForm * p , _gradientXForm * g , w , d ) ); }

		protected:
			// A reference to the underlying stream
			OutputLevelSetVertexStream< Real , Dim , Data > &_stream;

			// The affine transformation to be applied to the positions
			XForm< Real , Dim+1 > _positionXForm;

			// The linear transformation to be applied to the normals
			XForm< Real , Dim > _gradientXForm;
		};

		///////////////////////////////////////
		// Transformed Indexed Vertex Stream //
		///////////////////////////////////////
		template< typename Real , unsigned int Dim >
		struct TransformedOutputIndexedLevelSetVertexStream< Real , Dim > : public OutputIndexedLevelSetVertexStream< Real , Dim >
		{
			// A constructor initialized with the transformation to be applied to the samples, and a sample stream
			TransformedOutputIndexedLevelSetVertexStream( XForm< Real , Dim+1 > xForm , OutputIndexedLevelSetVertexStream< Real , Dim > &stream ) : _stream(stream) , _positionXForm(xForm)
			{
				_gradientXForm = XForm< Real , Dim > ( xForm ).inverse().transpose() * (Real)pow( xForm.determinant() , 1./Dim );
			}

			// Need to write the union to ensure that the counter gets set
			void base_write( node_index_type idx , Position< Real , Dim > p , Normal< Real , Dim > g , Real w )
			{
				_stream.write( std::pair< node_index_type , LevelSetVertex< Real , Dim > >( idx , LevelSetVertex< Real , Dim > (_positionXForm * p , _gradientXForm * g , w ) ) );
			}
			void base_write( unsigned int thread , node_index_type idx , Position< Real , Dim > p , Normal< Real , Dim > g , Real w )
			{
				_stream.write( thread , std::pair< node_index_type , LevelSetVertex< Real , Dim > >( idx , LevelSetVertex< Real , Dim > (_positionXForm * p , _gradientXForm * g , w ) ) );
			}

		protected:
			// A reference to the underlying stream
			OutputIndexedLevelSetVertexStream< Real , Dim > &_stream;

			// The affine transformation to be applied to the positions
			XForm< Real , Dim+1 > _positionXForm;

			// The linear transformation to be applied to the normals
			XForm< Real , Dim > _gradientXForm;
		};


		///////////////////////////////////////////////
		// Transformed Indexed Vertex w/ Data Stream //
		///////////////////////////////////////////////
		template< typename Real , unsigned int Dim , typename Data >
		struct TransformedOutputIndexedLevelSetVertexStream< Real , Dim , Data > : public OutputIndexedLevelSetVertexStream< Real , Dim , Data >
		{
			// A constructor initialized with the transformation to be applied to the samples, and a sample stream
			TransformedOutputIndexedLevelSetVertexStream( XForm< Real , Dim+1 > xForm , OutputIndexedLevelSetVertexStream< Real , Dim , Data > &stream ) : _stream(stream) , _positionXForm(xForm)
			{
				_gradientXForm = XForm< Real , Dim > ( xForm ).inverse().transpose() * (Real)pow( xForm.determinant() , 1./Dim );
			}

			void base_write( node_index_type idx , Position< Real , Dim > p , Normal< Real , Dim > g , Real w , Data d )
			{
				_stream.write( std::pair< node_index_type , LevelSetVertex< Real , Dim , Data > >( idx , LevelSetVertex< Real , Dim , Data >( _positionXForm * p , _gradientXForm * g , w , d ) ) );
			}
			void base_write( unsigned int thread , node_index_type idx , Position< Real , Dim > p , Normal< Real , Dim > g , Real w , Data d )
			{
				_stream.write( thread , std::pair< node_index_type , LevelSetVertex< Real , Dim , Data > >( idx , LevelSetVertex< Real , Dim , Data >( _positionXForm * p , _gradientXForm * g , w , d ) ) );
			}

		protected:
			// A reference to the underlying stream
			OutputIndexedLevelSetVertexStream< Real , Dim , Data > &_stream;

			// The affine transformation to be applied to the positions
			XForm< Real , Dim+1 > _positionXForm;

			// The linear transformation to be applied to the normals
			XForm< Real , Dim > _gradientXForm;
		};

		///////////////////////////////////////////
		// A wrapper class to write out vertices //
		///////////////////////////////////////////
		template< typename Vertex , typename Real , unsigned int Dim >
		struct OutputLevelSetVertexStreamWrapper< Vertex , Real , Dim > : public OutputLevelSetVertexStream< Real , Dim >
		{
			virtual Vertex toOutputVertex( const LevelSetVertex< Real , Dim > &in ) = 0;

			OutputLevelSetVertexStreamWrapper( OutputDataStream< Vertex > &stream ) : _stream(stream){}

			void base_write( Position< Real , Dim > p , Normal< Real , Dim > g , Real w )
			{
				LevelSetVertex< Real , Dim > in;
				in.template get<0>() = p;
				in.template get<1>() = g;
				in.template get<2>() = w;

				Vertex out = toOutputVertex( in );
				_stream.write( out );
			}
			void base_write( unsigned int thread , Position< Real , Dim > p , Normal< Real , Dim > g , Real w )
			{
				LevelSetVertex< Real , Dim > in;
				in.template get<0>() = p;
				in.template get<1>() = g;
				in.template get<2>() = w;

				Vertex out = toOutputVertex( in );
				_stream.write( thread , out );
			}
		protected:
			OutputDataStream< Vertex > &_stream;
		};

		template< typename Vertex , typename Real , unsigned int Dim , typename Data >
		struct OutputLevelSetVertexStreamWrapper< Vertex , Real , Dim , Data > : public OutputLevelSetVertexStream< Real , Dim , Data >
		{
			virtual Vertex toOutputVertex( const LevelSetVertex< Real , Dim , Data > &in ) = 0;

			OutputLevelSetVertexStreamWrapper( OutputDataStream< Vertex > &stream ) : _stream(stream){}

			void base_write( Position< Real , Dim > p , Normal< Real , Dim > g , Real w , Data d )
			{
				LevelSetVertex< Real , Dim , Data > in;
				in.template get<0>() = p;
				in.template get<1>() = g;
				in.template get<2>() = w;
				in.template get<3>() = d;

				Vertex out = toOutputVertex( in );
				_stream.write( out );
			}
			void base_write( unsigned int thread , Position< Real , Dim > p , Normal< Real , Dim > g , Real w , Data d )
			{
				LevelSetVertex< Real , Dim , Data > in;
				in.template get<0>() = p;
				in.template get<1>() = g;
				in.template get<2>() = w;
				in.template get<3>() = d;

				Vertex out = toOutputVertex( in );
				_stream.write( thread , out );
			}
		protected:
			OutputDataStream< Vertex > &_stream;
		};

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// A wrapper class to write out indexed vertices: OutputDataStream< IndexedVertex > -> OutputIndexedLevelSetVertexStream< Real , Dim > //
		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		template< typename Vertex , typename Real , unsigned int Dim >
		struct OutputIndexedLevelSetVertexStreamWrapper< Vertex , Real , Dim > : public OutputIndexedLevelSetVertexStream< Real , Dim >
		{
			virtual Vertex toOutputVertex( const LevelSetVertex< Real , Dim > &in ) = 0;

			OutputIndexedLevelSetVertexStreamWrapper( OutputDataStream< std::pair< node_index_type , Vertex > > &stream ) : _stream(stream){}

			void base_write( node_index_type idx , Position< Real , Dim > p , Normal< Real , Dim > g , Real w )
			{
				LevelSetVertex< Real , Dim > in;
				in.template get<0>() = p;
				in.template get<1>() = g;
				in.template get<2>() = w;

				std::pair< node_index_type , Vertex > out( idx , toOutputVertex(in) );
				_stream.write( out );
			}
			void base_write( unsigned int thread , node_index_type idx , Position< Real , Dim > p , Normal< Real , Dim > g , Real w )
			{
				LevelSetVertex< Real , Dim > in;
				in.template get<0>() = p;
				in.template get<1>() = g;
				in.template get<2>() = w;

				std::pair< node_index_type , Vertex > out( idx , toOutputVertex(in) );
				_stream.write( thread , out );
			}
		protected:
			OutputDataStream< std::pair< node_index_type , Vertex > > &_stream;
		};

		template< typename Vertex , typename Real , unsigned int Dim , typename Data >
		struct OutputIndexedLevelSetVertexStreamWrapper< Vertex , Real , Dim , Data > : public OutputIndexedLevelSetVertexStream< Real , Dim , Data >
		{
			virtual Vertex toOutputVertex( const LevelSetVertex< Real , Dim , Data > &in ) = 0;

			OutputIndexedLevelSetVertexStreamWrapper( OutputDataStream< std::pair< node_index_type , Vertex > > &stream ) : _stream(stream){}

			void base_write( node_index_type idx , Position< Real , Dim > p , Normal< Real , Dim > g , Real w , Data d )
			{
				LevelSetVertex< Real , Dim , Data > in;
				in.template get<0>() = p;
				in.template get<1>() = g;
				in.template get<2>() = w;
				in.template get<3>() = d;

				std::pair< node_index_type , Vertex > out( idx , toOutputVertex(in) );
				_stream.write( out );
			}
			void base_write( unsigned int thread , node_index_type idx , Position< Real , Dim > p , Normal< Real , Dim > g , Real w , Data d )
			{
				LevelSetVertex< Real , Dim , Data > in;
				in.template get<0>() = p;
				in.template get<1>() = g;
				in.template get<2>() = w;
				in.template get<3>() = d;

				std::pair< node_index_type , Vertex > out( idx , toOutputVertex(in) );
				_stream.write( thread , out );
			}
		protected:
			OutputDataStream< std::pair< node_index_type , Vertex > > &_stream;
		};
		//////////////////////////////////
		// File-backed streaming memory //
		//////////////////////////////////
		class FileBackedReadWriteStream
		{
		public:
			struct FileDescription
			{
				FILE *fp;

				FileDescription( FILE *fp ) : fp(fp) , _closeFile(false)
				{
					if( !this->fp )
					{
						this->fp = std::tmpfile();
						_closeFile = true;
						if( !this->fp ) ERROR_OUT( "Failed to open temporary file" );
					}
				}
				~FileDescription( void ){ if( _closeFile ) fclose(fp); }
			protected:
#ifdef SHOW_WARNINGS
#pragma message( "[WARNING] Probably can let the system handle closing the file" )
#endif // SHOW_WARNINGS
				bool _closeFile;
			};

			FileBackedReadWriteStream( FILE *fp ) : _fd(fp) {}
			bool write( ConstPointer(char) data , size_t size ){ return fwrite( data , sizeof(char) , size , _fd.fp )==size; }
			bool read( Pointer(char) data , size_t size ){ return fread( data , sizeof(char) , size , _fd.fp )==size; }
			void reset( void ){ fseek( _fd.fp , 0 , SEEK_SET ); }
		protected:
			FileDescription _fd;
		};

		//////////////////////////////////////////////////////////////////////////////
		// Output and the input face stream, backed either by a file or by a vector //
		//////////////////////////////////////////////////////////////////////////////
		// [WARNING] These assume that the stream starts as write-only and after the reset method is invoked, the stream becomes read-only.

		template< unsigned int FaceDim >
		struct OutputInputFaceStream : public OutputFaceStream< FaceDim > , public InputFaceStream< FaceDim >
		{
			// The streams for communicating the information
			InputFaceStream < FaceDim > * inStream;
			OutputFaceStream< FaceDim > *outStream;

			void reset( void ){ inStream->reset(); }
			bool base_read( Face< FaceDim > &f ){ return inStream->read(f); }
			bool base_read( unsigned int t , Face< FaceDim > &f ){ return inStream->read(t,f); }
			void base_write( const Face< FaceDim > &f ){ outStream->write(f); }
			void base_write( unsigned int t , const Face< FaceDim > &f ){ outStream->write(t,f); }

			OutputInputFaceStream( bool inCore , bool multi )
			{
				size_t sz = std::thread::hardware_concurrency();

				_backingVector = NULL;
				_backingVectors.resize( sz , NULL );

				_backingFile = NULL;
				_backingFiles.resize( sz , NULL );

				_inStreams.resize( sz , NULL );
				_outStreams.resize( sz , NULL );

				if( inCore )
				{
					if( multi )
					{
						for( unsigned int i=0 ; i<sz ; i++ )
						{
							_backingVectors[i] = new std::vector< Face< FaceDim > >();
							_inStreams[i] = new VectorBackedInputDataStream< Face< FaceDim > >( *_backingVectors[i] );
							_outStreams[i] = new VectorBackedOutputDataStream< Face< FaceDim > >( *_backingVectors[i] );
						}
						inStream = new MultiInputDataStream< Face< FaceDim > >( _inStreams );
						outStream = new MultiOutputDataStream< Face< FaceDim > >( _outStreams );
					}
					else
					{
						_backingVector = new std::vector< Face< FaceDim > >();
						inStream = new VectorBackedInputDataStream< Face< FaceDim > >( *_backingVector );
						outStream = new VectorBackedOutputDataStream< Face< FaceDim > >( *_backingVector );
					}
				}
				else
				{
					if( multi )
					{
						for( unsigned int i=0 ; i<sz ; i++ )
						{
							_backingFiles[i] = new FileBackedReadWriteStream::FileDescription( NULL );
							_inStreams[i] = new FileBackedInputDataStream< Face< FaceDim > >( _backingFiles[i]->fp );
							_outStreams[i] = new FileBackedOutputDataStream< Face< FaceDim > >( _backingFiles[i]->fp );
						}
						inStream = new MultiInputDataStream< Face< FaceDim > >( _inStreams );
						outStream = new MultiOutputDataStream< Face< FaceDim > >( _outStreams );
					}
					else
					{
						_backingFile = new FileBackedReadWriteStream::FileDescription( NULL );
						inStream = new FileBackedInputDataStream< Face< FaceDim > >( _backingFile->fp );
						outStream = new FileBackedOutputDataStream< Face< FaceDim > >( _backingFile->fp );
					}
				}
			}

			~OutputInputFaceStream( void )
			{
				size_t sz = std::thread::hardware_concurrency();

				delete _backingVector;
				delete _backingFile;

				for( unsigned int i=0 ; i<sz ; i++ )
				{
					delete _backingVectors[i];
					delete _backingFiles[i];
					delete  _inStreams[i];
					delete _outStreams[i];
				}

				delete  inStream;
				delete outStream;
			}
		protected:
			std::vector< Face< FaceDim > > *_backingVector;
			FileBackedReadWriteStream::FileDescription *_backingFile;
			std::vector< std::vector< Face< FaceDim > > * > _backingVectors;
			std::vector< FileBackedReadWriteStream::FileDescription * > _backingFiles;
			std::vector<  InputDataStream< Face< FaceDim > > * >  _inStreams;
			std::vector< OutputDataStream< Face< FaceDim > > * > _outStreams;
		};

		template< typename Factory , bool Parallel >
		struct OutputInputFactoryTypeStream :
			public std::conditional_t< Parallel , MultiOutputDataStream< std::pair< node_index_type , typename Factory::VertexType > > , OutputDataStream< typename Factory::VertexType > > ,
			public std::conditional_t< Parallel ,  MultiInputDataStream< std::pair< node_index_type , typename Factory::VertexType > > ,  InputDataStream< typename Factory::VertexType > >
		{
			using Vertex = std::conditional_t< Parallel , std::pair< node_index_type , typename Factory::VertexType > , typename Factory::VertexType >;
			// The streams for communicating the information
			using  InputStreamType = std::conditional_t< Parallel ,  MultiInputDataStream< Vertex > ,  InputDataStream< Vertex > >;
			using OutputStreamType = std::conditional_t< Parallel , MultiOutputDataStream< Vertex > , OutputDataStream< Vertex > >;
			InputStreamType * inStream;
			OutputStreamType *outStream;

			void reset( void ){ inStream->reset(); }
			void base_write(                       const Vertex &v ){ outStream->write(          v ); }
			void base_write( unsigned int thread , const Vertex &v ){ outStream->write( thread , v ); }
			bool base_read(                       Vertex &v ){ return inStream->read(          v ); }
			bool base_read( unsigned int thread , Vertex &v ){ return inStream->read( thread , v ); }

			OutputInputFactoryTypeStream( Factory &factory , bool inCore )
			{
				size_t sz = std::thread::hardware_concurrency();

				_backingVector = NULL;
				_backingVectors.resize( sz , NULL );

				_backingFile = NULL;
				_backingFiles.resize( sz , NULL );

				_inStreams.resize( sz , NULL );
				_outStreams.resize( sz , NULL );

				if( inCore )
				{
					if constexpr( Parallel )
					{
						for( unsigned int i=0 ; i<sz ; i++ )
						{
							_backingVectors[i] = new std::vector< Vertex >();
							_inStreams[i] = new VectorBackedInputDataStream< Vertex >( *_backingVectors[i] );
							_outStreams[i] = new VectorBackedOutputDataStream< Vertex >( *_backingVectors[i] );
						}
						inStream = new MultiInputDataStream< Vertex >( _inStreams );
						outStream = new MultiOutputDataStream< Vertex >( _outStreams );
					}
					else
					{
						_backingVector = new std::vector< Vertex >();

						inStream = new VectorBackedInputDataStream< Vertex >( *_backingVector );
						outStream = new VectorBackedOutputDataStream< Vertex >( *_backingVector );
					}
				}
				else
				{
					if constexpr( Parallel )
					{
						for( unsigned int i=0 ; i<sz ; i++ )
						{
							_backingFiles[i] = new FileBackedReadWriteStream::FileDescription( NULL );
							_inStreams[i] = new FileBackedInputFactoryTypeStream< Factory , Parallel , node_index_type >( _backingFiles[i]->fp , factory );
							_outStreams[i] = new FileBackedOutputFactoryTypeStream< Factory , Parallel , node_index_type >( _backingFiles[i]->fp , factory  );
						}
						inStream = new MultiInputDataStream< Vertex >( _inStreams );
						outStream = new MultiOutputDataStream< Vertex >( _outStreams );
					}
					else
					{
						_backingFile = new FileBackedReadWriteStream::FileDescription( NULL );
						inStream = new FileBackedInputFactoryTypeStream< Factory , Parallel , node_index_type >( _backingFile->fp , factory );
						outStream = new FileBackedOutputFactoryTypeStream< Factory , Parallel , node_index_type >( _backingFile->fp , factory );
					}
				}
				if constexpr( Parallel )
				{
					MultiOutputDataStream< std::pair< node_index_type , typename Factory::VertexType > >::_init( _outStreams );
					MultiInputDataStream< std::pair< node_index_type , typename Factory::VertexType > >::_init( _inStreams );
				}
			}

			~OutputInputFactoryTypeStream( void )
			{
				size_t sz = std::thread::hardware_concurrency();

				delete _backingVector;
				delete _backingFile;

				for( unsigned int i=0 ; i<sz ; i++ )
				{
					delete _backingVectors[i];
					delete _backingFiles[i];
					delete  _inStreams[i];
					delete _outStreams[i];
				}

				delete  inStream;
				delete outStream;
			}
		protected:
			std::vector< Vertex > *_backingVector;
			FileBackedReadWriteStream::FileDescription *_backingFile;
			std::vector< std::vector< Vertex > * >_backingVectors;
			std::vector< FileBackedReadWriteStream::FileDescription * > _backingFiles;
			std::vector<  InputDataStream< Vertex > * >  _inStreams;
			std::vector< OutputDataStream< Vertex > * > _outStreams;
		};

		template< typename Real , unsigned int Dim , bool HasGradients , bool HasDensity >
		struct OutputLevelSetVertexInfo< Real , Dim , HasGradients , HasDensity >
		{
			using Factory =
				typename std::conditional
				<
				HasGradients ,
				typename std::conditional
				<
				HasDensity ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::NormalFactory< Real , Dim > , VertexFactory::ValueFactory< Real > > ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::NormalFactory< Real , Dim > >
				>::type ,
				typename std::conditional
				<
				HasDensity ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::ValueFactory< Real > > ,
				VertexFactory::PositionFactory< Real , Dim >
				>::type
				>::type;
			using Vertex = typename Factory::VertexType;

			static Factory GetFactory( void ){ return Factory(); }

			struct StreamWrapper : public Reconstructor::OutputLevelSetVertexStreamWrapper< Vertex , Real , Dim >
			{
				StreamWrapper( OutputDataStream< Vertex > &stream ) :
					Reconstructor::OutputLevelSetVertexStreamWrapper< Vertex , Real , Dim >( stream ){}

				Vertex toOutputVertex( const Reconstructor::LevelSetVertex< Real , Dim > &in )
				{
					Vertex out;
					if constexpr( HasGradients || HasDensity )
					{
						out.template get<0>() = in.template get<0>();
						if constexpr( HasGradients )
						{
							out.template get<1>() = in.template get<1>();
							if constexpr( HasDensity ) out.template get<2>() = in.template get<2>();
						}
						else
						{
							if constexpr( HasDensity ) out.template get<1>() = in.template get<2>();
						}
					}
					else out = in.template get<0>();
					return out;
				}
			};
		};

		template< typename Real , unsigned int Dim , bool HasGradients , bool HasDensity , typename AuxDataFactory >
		struct OutputLevelSetVertexInfo< Real , Dim , HasGradients , HasDensity , AuxDataFactory >
		{
			using Factory =
				typename std::conditional
				<
				HasGradients ,
				typename std::conditional
				<
				HasDensity ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::NormalFactory< Real , Dim > , VertexFactory::ValueFactory< Real > , AuxDataFactory > ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::NormalFactory< Real , Dim > , AuxDataFactory >
				>::type ,
				typename std::conditional
				<
				HasDensity ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::ValueFactory< Real > , AuxDataFactory > ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , AuxDataFactory >
				>::type
				>::type;
			using AuxData = typename AuxDataFactory::VertexType;

			using _Vertex = VectorTypeUnion< Real , Point< Real , Dim > , Point< Real , Dim > , Real , typename AuxDataFactory::VertexType >;
			using Vertex = typename Factory::VertexType;

			static Factory GetFactory( AuxDataFactory auxDataFactory )
			{
				if constexpr( HasGradients )
				{
					if constexpr( HasDensity ) return Factory( VertexFactory::PositionFactory< Real , Dim >() , VertexFactory::NormalFactory< Real , Dim >() , VertexFactory::ValueFactory< Real >() , auxDataFactory );
					else                       return Factory( VertexFactory::PositionFactory< Real , Dim >() , VertexFactory::NormalFactory< Real , Dim >() ,                                         auxDataFactory );
				}
				else
				{
					if constexpr( HasDensity ) return Factory( VertexFactory::PositionFactory< Real , Dim >() ,                                                VertexFactory::ValueFactory< Real >() , auxDataFactory );
					else                       return Factory( VertexFactory::PositionFactory< Real , Dim >() ,                                                                                        auxDataFactory );
				}
			}

			struct StreamWrapper : public Reconstructor::OutputLevelSetVertexStreamWrapper< Vertex , Real , Dim , AuxData >
			{
				StreamWrapper( OutputDataStream< Vertex > &stream ) :
					Reconstructor::OutputLevelSetVertexStreamWrapper< Vertex , Real , Dim , AuxData >( stream ){}

				Vertex toOutputVertex( const Reconstructor::LevelSetVertex< Real , Dim , AuxData > &in )
				{
					Vertex out;
					out.template get<0>() = in.template get<0>();
					if constexpr( HasGradients )
					{
						out.template get<1>() = in.template get<1>();
						if constexpr( HasDensity ) out.template get<2>() = in.template get<2>() , out.template get<3>() = in.template get<3>();
						else out.template get<2>() = in.template get<3>();
					}
					else
					{
						if constexpr( HasDensity ) out.template get<1>() = in.template get<2>() , out.template get<2>() = in.template get<3>();
						else out.template get<1>() = in.template get<3>();
					}
					return out;
				}
			};
		};
		template< typename Real , unsigned int Dim , bool HasGradients , bool HasDensity >
		struct OutputIndexedLevelSetVertexInfo< Real , Dim , HasGradients , HasDensity >
		{
			using Factory =
				typename std::conditional
				<
				HasGradients ,
				typename std::conditional
				<
				HasDensity ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::NormalFactory< Real , Dim > , VertexFactory::ValueFactory< Real > > ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::NormalFactory< Real , Dim > >
				>::type ,
				typename std::conditional
				<
				HasDensity ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::ValueFactory< Real > > ,
				VertexFactory::PositionFactory< Real , Dim >
				>::type
				>::type;
			using Vertex = typename Factory::VertexType;
			using IndexedVertex = std::pair< node_index_type , Vertex >;

			static Factory GetFactory( void ){ return Factory(); }

			struct StreamWrapper : public Reconstructor::OutputIndexedLevelSetVertexStreamWrapper< Vertex , Real , Dim >
			{
				StreamWrapper( OutputDataStream< IndexedVertex > &stream ) :
					Reconstructor::OutputIndexedLevelSetVertexStreamWrapper< Vertex , Real , Dim >( stream ){}

				Vertex toOutputVertex( const Reconstructor::LevelSetVertex< Real , Dim > &in )
				{
					Vertex out;
					if constexpr( HasGradients || HasDensity )
					{
						out.template get<0>() = in.template get<0>();
						if constexpr( HasGradients )
						{
							out.template get<1>() = in.template get<1>();
							if constexpr( HasDensity ) out.template get<2>() = in.template get<2>();
						}
						else
						{
							if constexpr( HasDensity ) out.template get<1>() = in.template get<2>();
						}
					}
					else out = in.template get<0>();
					return out;
				}
			};
		};

		template< typename Real , unsigned int Dim , bool HasGradients , bool HasDensity , typename AuxDataFactory >
		struct OutputIndexedLevelSetVertexInfo< Real , Dim , HasGradients , HasDensity , AuxDataFactory >
		{
			using Factory =
				typename std::conditional
				<
				HasGradients ,
				typename std::conditional
				<
				HasDensity ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::NormalFactory< Real , Dim > , VertexFactory::ValueFactory< Real > , AuxDataFactory > ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::NormalFactory< Real , Dim > , AuxDataFactory >
				>::type ,
				typename std::conditional
				<
				HasDensity ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , VertexFactory::ValueFactory< Real > , AuxDataFactory > ,
				VertexFactory::Factory< Real , VertexFactory::PositionFactory< Real , Dim > , AuxDataFactory >
				>::type
				>::type;
			using AuxData = typename AuxDataFactory::VertexType;

			using _Vertex = VectorTypeUnion< Real , Point< Real , Dim > , Point< Real , Dim > , Real , typename AuxDataFactory::VertexType >;
			using Vertex = typename Factory::VertexType;
			using IndexedVertex = std::pair< node_index_type , Vertex >;

			static Factory GetFactory( AuxDataFactory auxDataFactory )
			{
				if constexpr( HasGradients )
				{
					if constexpr( HasDensity ) return Factory( VertexFactory::PositionFactory< Real , Dim >() , VertexFactory::NormalFactory< Real , Dim >() , VertexFactory::ValueFactory< Real >() , auxDataFactory );
					else                       return Factory( VertexFactory::PositionFactory< Real , Dim >() , VertexFactory::NormalFactory< Real , Dim >() ,                                         auxDataFactory );
				}
				else
				{
					if constexpr( HasDensity ) return Factory( VertexFactory::PositionFactory< Real , Dim >() ,                                                VertexFactory::ValueFactory< Real >() , auxDataFactory );
					else                       return Factory( VertexFactory::PositionFactory< Real , Dim >() ,                                                                                        auxDataFactory );
				}
			}

			struct StreamWrapper : public Reconstructor::OutputIndexedLevelSetVertexStreamWrapper< Vertex , Real , Dim , AuxData >
			{
				StreamWrapper( OutputDataStream< IndexedVertex > &stream ) :
					Reconstructor::OutputIndexedLevelSetVertexStreamWrapper< Vertex , Real , Dim , AuxData >( stream ){}

				Vertex toOutputVertex( const Reconstructor::LevelSetVertex< Real , Dim , AuxData > &in )
				{
					Vertex out;
					out.template get<0>() = in.template get<0>();
					if constexpr( HasGradients )
					{
						out.template get<1>() = in.template get<1>();
						if constexpr( HasDensity ) out.template get<2>() = in.template get<2>() , out.template get<3>() = in.template get<3>();
						else out.template get<2>() = in.template get<3>();
					}
					else
					{
						if constexpr( HasDensity ) out.template get<1>() = in.template get<2>() , out.template get<2>() = in.template get<3>();
						else out.template get<1>() = in.template get<3>();
					}
					return out;
				}
			};
		};
	}
}

#endif // RECONSTRUCTORS_STREAMS_INCLUDED