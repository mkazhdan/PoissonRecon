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

#undef ARRAY_DEBUG

#include <tinyxml2.h>
#include <iostream>
#include <string>
#include "boost/format.hpp"

#include "Mesh/PoissonRecon/SurfaceTrimming.h"

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

void usage()
{
  std::cerr << " Usage: ./exec <path/to/xml> <start_frame_number> "
               "<end_frame_number> <calculate_normal>"
            << std::endl;
}

void ReadXML(tinyxml2::XMLDocument &doc_xml,
             std::string &input_path,
             std::string &output_path,
             std::string &trim,
             std::string &smooth,
             std::string &a_ratio,
             std::string &polygon_mesh)
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

  tinyxml2::XMLHandle trim_xml = input_xml.FirstChildElement("TRIM");

  if (trim_xml.ToElement())
  {
    trim = trim_xml.ToElement()->GetText();

    std::cout << "read trim" << std::endl;
  }

  tinyxml2::XMLHandle smooth_xml = input_xml.FirstChildElement("SMOOTH");

  if (smooth_xml.ToElement())
  {
    smooth = smooth_xml.ToElement()->GetText();

    std::cout << "read smooth" << std::endl;
  }

  tinyxml2::XMLHandle a_ratio_xml = input_xml.FirstChildElement("A_RATIO");

  if (a_ratio_xml.ToElement())
  {
    a_ratio = a_ratio_xml.ToElement()->GetText();

    std::cout << "read a_ratio" << std::endl;
  }

  tinyxml2::XMLHandle polygon_mesh_xml =
      input_xml.FirstChildElement("POLYGON_MESH");

  if (polygon_mesh_xml.ToElement())
  {
    polygon_mesh = polygon_mesh_xml.ToElement()->GetText();

    std::cout << "read polygon_mesh" << std::endl;
  }

  tinyxml2::XMLHandle output_xml =
      xml.FirstChildElement("CONFIG").FirstChildElement("OUTPUT");

  tinyxml2::XMLHandle output_path_xml =
      output_xml.FirstChildElement("OUTPUT_PATH");
  if (output_path_xml.ToElement())
  {
    output_path = output_path_xml.ToElement()->GetText();
    std::cout << "read output path" << std::endl;
  }
}

int main(int argc, char **argv)
{
  if (argc != 4)
  {
    usage();
    return -1;
  }

  std::string input_path;
  std::string output_path;
  std::string trim;
  std::string smooth;
  std::string a_ratio;
  std::string polygon_mesh;

  tinyxml2::XMLDocument doc_xml;
  if (doc_xml.LoadFile(argv[1]))
  {
    std::cerr << "failed to load xml file" << std::endl;
    return 0;
  }
  else
  {
    ReadXML(
        doc_xml, input_path, output_path, trim, smooth, a_ratio, polygon_mesh);
  }
  const std::size_t start_frm = std::atoi(argv[2]);
  const std::size_t end_frm   = std::atoi(argv[3]);
  for (size_t k = start_frm; k <= end_frm; k++)
  {
    std::string temp_in  = ParsePath(input_path, k);
    std::string temp_out = ParsePath(output_path, k);

    std::cout << temp_in << std::endl;
    std::cout << temp_out << std::endl;
    char *options[] = {const_cast<char *>("--in"),
                       const_cast<char *>(temp_in.c_str()),
                       const_cast<char *>("--out"),
                       const_cast<char *>(temp_out.c_str()),
                       const_cast<char *>("--trim"),
                       const_cast<char *>(trim.c_str()),
                       const_cast<char *>("--smooth"),
                       const_cast<char *>(smooth.c_str()),
                       const_cast<char *>("--aRatio"),
                       const_cast<char *>(a_ratio.c_str()),
                       const_cast<char *>(polygon_mesh.c_str())};

    int nOptions = sizeof(options) / sizeof(char *);

    SurfaceTrimming trim;
    trim.trim(nOptions, options);
  }
}
