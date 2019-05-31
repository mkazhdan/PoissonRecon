/*
Copyright (c) 2006, Michael Kazhdan and Matthew Bolitho
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list
of conditions and the following disclaimer. Redistributions in binary form must
reproduce the above copyright notice, this list of conditions and the following
disclaimer in the documentation and/or other materials provided with the
distribution.

Neither the name of the Johns Hopkins University nor the names of its
contributors may be used to endorse or promote products derived from this
software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <tinyxml2.h>
#include <iostream>
#include <string>

#include "Mesh/PoissonRecon/PoissonRecon.h"

template <class... Int>
std::string ParsePath(const std::string &path, const Int... index)
{
  if (path.find('%') != std::string::npos)
  {
    char buffer[4096];
    return std::string(
        buffer, std::snprintf(buffer, sizeof(buffer), path.c_str(), index...));
  }

  return path;
}

void usage()
{
  std::cerr << " Usage: ./exec <path/to/xml> <start_frame_number> "
               "<end_frame_number> "
            << std::endl;
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
             std::string &iters,
             std::string &primal_voxel,
             std::string &colors,
             std::string &data,
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

  tinyxml2::XMLHandle iters_xml = input_xml.FirstChildElement("ITERS");

  if (iters_xml.ToElement())
  {
    iters = iters_xml.ToElement()->GetText();

    std::cout << "read iters" << std::endl;
  }

  tinyxml2::XMLHandle primal_voxel_xml =
      input_xml.FirstChildElement("PRIMAL_VOXEL");

  if (primal_voxel_xml.ToElement())
  {
    primal_voxel = primal_voxel_xml.ToElement()->GetText();

    std::cout << "primal_voxel_read" << std::endl;
  }

  tinyxml2::XMLHandle colors_xml = input_xml.FirstChildElement("COLORS");

  if (colors_xml.ToElement())
  {
    colors = colors_xml.ToElement()->GetText();

    std::cout << "read colors" << std::endl;
  }

  tinyxml2::XMLHandle data_xml = input_xml.FirstChildElement("DATA");

  if (data_xml.ToElement())
  {
    data = data_xml.ToElement()->GetText();

    std::cout << "read data" << std::endl;
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

int main(int argc, char *argv[])
{
  if (argc != 4)
  {
    usage();
    return -1;
  }

  std::string input_path;
  std::string output_path;
  std::string degree;
  std::string b_type;
  std::string depth;
  std::string scale;
  std::string samples_per_node;
  std::string point_weight;
  std::string iters;
  std::string primal_voxel;
  std::string colors;
  std::string data;
  std::string density;
  std::string linear_fit;
  std::string polygon_mesh;
  std::string threads;
  std::string verbose;
  std::string temp = "/home/hypevr/Poisson/PoissonRecon/";

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
            iters,
            primal_voxel,
            colors,
            data,
            density,
            linear_fit,
            polygon_mesh,
            threads,
            verbose);
  }
  const std::size_t start_frm = std::atoi(argv[2]);
  const std::size_t end_frm   = std::atoi(argv[3]);
  for (size_t k = start_frm; k <= end_frm; k++)
  {
    std::string temp_in  = ParsePath(input_path, k);
    std::string temp_out = ParsePath(output_path, k);

    std::cout << temp_in << std::endl;
    std::cout << temp_out << std::endl;

    char *options[] = {const_cast<char *>("PoissonRecon"),
                       const_cast<char *>("--in"),
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
                       const_cast<char *>("--iters"),
                       const_cast<char *>(iters.c_str()),
                       const_cast<char *>(primal_voxel.c_str()),
                       const_cast<char *>(colors.c_str()),
                       const_cast<char *>("--data"),
                       const_cast<char *>(data.c_str()),
                       const_cast<char *>(density.c_str()),
                       const_cast<char *>(linear_fit.c_str()),
                       const_cast<char *>(polygon_mesh.c_str()),
                       const_cast<char *>("--threads"),
                       const_cast<char *>(threads.c_str()),
                       const_cast<char *>("--tempDir"),
                       const_cast<char *>(temp.c_str()),
                       const_cast<char *>(verbose.c_str()),
                       const_cast<char *>("--ascii")};

    int nOptions = sizeof(options) / sizeof(char *);

    PoissonRecon poisson;
    poisson.reconstruct(nOptions, options);
  }
}
