/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of
conditions and the following disclaimer. Redistributions in binary form must
reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution.

Neither the name of the Johns Hopkins University nor the names of its
contributors
may be used to endorse or promote products derived from this software without
specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED
WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
SUCH
DAMAGE.
*/

#undef FAST_COMPILE
#undef ARRAY_DEBUG
#define BRUNO_LEVY_FIX
#define FOR_RELEASE

#include <tinyxml2.h>
#include <string>
#include "boost/format.hpp"

#include "Mesh/PoissonRecon/PoissonRecon.h"

std::string ParsePath(std::string path, const int index)
{
  if (path.find('%') != std::string::npos)
  {
    return boost::str(boost::format(path) % index);
  }
  else
  {
    return path;
  }
}

void ReadXML(tinyxml2::XMLDocument &doc_xml,
             std::string &input_path,
             std::string &output_path,
             std::string &degree,
             std::string &b_type,
             std::string &depth,
             std::string &scale,
             std::string &samples_per_node,
             std::string &point_weight,
             std::string &confidence,
             std::string &n_weights,
             std::string &iters,
             std::string &cg_depth,
             std::string &full_depth,
             std::string &voxel_depth,
             std::string &primal_voxel,
             std::string &color,
             std::string &density,
             std::string &linear_fit,
             std::string &polygon_mesh,
             std::string &threads,
             std::string &verbose)
{
  tinyxml2::XMLHandle xml(&doc_xml);

  tinyxml2::XMLHandle input_xml =
      xml.FirstChildElement("CONFIG").FirstChildElement("INPUT");

  tinyxml2::XMLHandle input_path_xml =
      input_xml.FirstChildElement("INPUT_PATH");
  if (input_path_xml.ToElement())
  {
    input_path = input_path_xml.ToElement()->GetText();
    std::cout << "read input path" << std::endl;
  }
  else
  {
    std::cerr << "NO INPUT SPECIFIED" << std::endl;
  }

  tinyxml2::XMLHandle degree_xml = input_xml.FirstChildElement("DEGREE");

  if (degree_xml.ToElement())
  {
    degree = degree_xml.ToElement()->GetText();

    std::cout << "read degree" << std::endl;
  }

  tinyxml2::XMLHandle b_type_xml = input_xml.FirstChildElement("B_TYPE");

  if (b_type_xml.ToElement())
  {
    b_type = b_type_xml.ToElement()->GetText();

    std::cout << "read b_type" << std::endl;
  }

  tinyxml2::XMLHandle depth_xml = input_xml.FirstChildElement("DEPTH");

  if (depth_xml.ToElement())
  {
    depth = depth_xml.ToElement()->GetText();

    std::cout << "read depth" << std::endl;
  }

  tinyxml2::XMLHandle scale_xml = input_xml.FirstChildElement("SCALE");

  if (scale_xml.ToElement())
  {
    scale = scale_xml.ToElement()->GetText();

    std::cout << "read scale" << std::endl;
  }

  tinyxml2::XMLHandle samples_per_node_xml =
      input_xml.FirstChildElement("SAMPLES_PER_NODE");

  if (samples_per_node_xml.ToElement())
  {
    samples_per_node = samples_per_node_xml.ToElement()->GetText();

    std::cout << "read samples_per_node" << std::endl;
  }

  tinyxml2::XMLHandle point_weight_xml =
      input_xml.FirstChildElement("POINT_WEIGHT");

  if (point_weight_xml.ToElement())
  {
    point_weight = point_weight_xml.ToElement()->GetText();

    std::cout << "read point_weight" << std::endl;
  }

  tinyxml2::XMLHandle confidence_xml =
      input_xml.FirstChildElement("CONFIDENCE");

  if (confidence_xml.ToElement())
  {
    confidence = confidence_xml.ToElement()->GetText();
    std::cout << "read confidence" << std::endl;
  }

  tinyxml2::XMLHandle n_weights_xml = input_xml.FirstChildElement("N_WEIGHTS");

  if (n_weights_xml.ToElement())
  {
    n_weights = n_weights_xml.ToElement()->GetText();
    std::cout << "read n_weights" << std::endl;
  }

  tinyxml2::XMLHandle iters_xml = input_xml.FirstChildElement("ITERS");

  if (iters_xml.ToElement())
  {
    iters = iters_xml.ToElement()->GetText();

    std::cout << "read iters" << std::endl;
  }

  tinyxml2::XMLHandle cg_depth_xml = input_xml.FirstChildElement("CG_DEPTH");

  if (cg_depth_xml.ToElement())
  {
    cg_depth = cg_depth_xml.ToElement()->GetText();

    std::cout << "read cg_depth" << std::endl;
  }

  tinyxml2::XMLHandle full_depth_xml =
      input_xml.FirstChildElement("FULL_DEPTH");

  if (full_depth_xml.ToElement())
  {
    full_depth = full_depth_xml.ToElement()->GetText();

    std::cout << "read full_depth" << std::endl;
  }

  tinyxml2::XMLHandle voxel_depth_xml =
      input_xml.FirstChildElement("VOXEL_DEPTH");

  if (voxel_depth_xml.ToElement())
  {
    voxel_depth = voxel_depth_xml.ToElement()->GetText();

    std::cout << "read voxel_depth" << std::endl;
  }

  tinyxml2::XMLHandle primal_voxel_xml =
      input_xml.FirstChildElement("PRIMAL_VOXEL");

  if (primal_voxel_xml.ToElement())
  {
    primal_voxel = primal_voxel_xml.ToElement()->GetText();

    std::cout << "primal_voxel_read" << std::endl;
  }

  tinyxml2::XMLHandle color_xml = input_xml.FirstChildElement("COLOR");

  if (color_xml.ToElement())
  {
    color = color_xml.ToElement()->GetText();

    std::cout << "read color" << std::endl;
  }

  tinyxml2::XMLHandle density_xml = input_xml.FirstChildElement("DENSITY");

  if (density_xml.ToElement())
  {
    density = density_xml.ToElement()->GetText();

    std::cout << "read density" << std::endl;
  }

  tinyxml2::XMLHandle linear_fit_xml =
      input_xml.FirstChildElement("LINEAR_FIT");

  if (linear_fit_xml.ToElement())
  {
    linear_fit = linear_fit_xml.ToElement()->GetText();

    std::cout << "read linear_fit" << std::endl;
  }

  tinyxml2::XMLHandle polygon_mesh_xml =
      input_xml.FirstChildElement("POLYGON_MESH");

  if (polygon_mesh_xml.ToElement())
  {
    polygon_mesh = polygon_mesh_xml.ToElement()->GetText();

    std::cout << "read polygon_mesh" << std::endl;
  }

  tinyxml2::XMLHandle threads_xml = input_xml.FirstChildElement("THREADS");

  if (threads_xml.ToElement())
  {
    threads = threads_xml.ToElement()->GetText();

    std::cout << "read threads" << std::endl;
  }

  tinyxml2::XMLHandle verbose_xml = input_xml.FirstChildElement("VERBOSE");

  if (verbose_xml.ToElement())
  {
    verbose = verbose_xml.ToElement()->GetText();

    std::cout << "will show verbose output" << std::endl;
  }

  tinyxml2::XMLHandle output_xml =
      xml.FirstChildElement("CONFIG").FirstChildElement("OUTPUT");

  tinyxml2::XMLHandle output_path_xml =
      output_xml.FirstChildElement("OUTPUT_PATH");

  if (output_path_xml.ToElement())
  {
    output_path = output_path_xml.ToElement()->GetText();

    std::cout << "read output_path" << std::endl;
  }

  else
  {
    std::cerr << "NO OUTPUT_PATH SPECIFIED" << std::endl;
  }
}

int main(int argc, char **argv)
{
  std::string input_path;
  std::string output_path;
  std::string degree;
  std::string b_type;
  std::string depth;
  std::string scale;
  std::string samples_per_node;
  std::string point_weight;
  std::string confidence;
  std::string n_weights;
  std::string iters;
  std::string cg_depth;
  std::string full_depth;
  std::string voxel_depth;
  std::string primal_voxel;
  std::string color;
  std::string density;
  std::string linear_fit;
  std::string polygon_mesh;
  std::string threads;
  std::string verbose;

  tinyxml2::XMLDocument doc_xml;
  if (doc_xml.LoadFile(argv[1]))
  {
    std::cerr << "failed to load xml file" << std::endl;
    return 0;
  }
  else
  {
    ReadXML(doc_xml,
            input_path,
            output_path,
            degree,
            b_type,
            depth,
            scale,
            samples_per_node,
            point_weight,
            confidence,
            n_weights,
            iters,
            cg_depth,
            full_depth,
            voxel_depth,
            primal_voxel,
            color,
            density,
            linear_fit,
            polygon_mesh,
            threads,
            verbose);
    std::cout << input_path << std::endl;
    std::cout << input_path.c_str() << std::endl;
    std::cout << ParsePath(input_path, 401) << std::endl;

    std::cout << ParsePath(input_path, 401).c_str() << std::endl;
    std::cout << const_cast<char *>(ParsePath(input_path, 401).c_str())
              << std::endl;
  }
  std::string temp_in  = ParsePath(input_path, 401);
  std::string temp_out = ParsePath(output_path, 401);

  char *options[] = {const_cast<char *>("--in"),
                     const_cast<char *>(temp_in.c_str()),
                     const_cast<char *>("--out"),
                     const_cast<char *>(temp_out.c_str()),
                     const_cast<char *>("--degree"),
                     const_cast<char *>(degree.c_str()),
                     const_cast<char *>("--bType"),
                     const_cast<char *>(b_type.c_str()),
                     const_cast<char *>("--depth"),
                     const_cast<char *>(depth.c_str()),
                     const_cast<char *>("--scale"),
                     const_cast<char *>(scale.c_str()),
                     const_cast<char *>("--samplesPerNode"),
                     const_cast<char *>(samples_per_node.c_str()),
                     const_cast<char *>("--pointWeight"),
                     const_cast<char *>(point_weight.c_str()),
                     const_cast<char *>(confidence.c_str()),
                     const_cast<char *>(n_weights.c_str()),
                     const_cast<char *>("--iters"),
                     const_cast<char *>(iters.c_str()),
                     const_cast<char *>("--cgDepth"),
                     const_cast<char *>(cg_depth.c_str()),
                     const_cast<char *>("--fullDepth"),
                     const_cast<char *>(full_depth.c_str()),
                     const_cast<char *>("--voxelDepth"),
                     const_cast<char *>(voxel_depth.c_str()),
                     const_cast<char *>(primal_voxel.c_str()),
                     const_cast<char *>("--color"),
                     const_cast<char *>(color.c_str()),
                     const_cast<char *>(density.c_str()),
                     const_cast<char *>(linear_fit.c_str()),
                     const_cast<char *>(polygon_mesh.c_str()),
                     const_cast<char *>("--threads"),
                     const_cast<char *>(threads.c_str()),
                     const_cast<char *>(verbose.c_str())};

  int nOptions = sizeof(options) / sizeof(char *);

  PoissonRecon poisson;

  poisson.reconstruct(nOptions, options);
}
