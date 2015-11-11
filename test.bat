\Research\PoissonRecon\PoissonRecon\Bin\x64\Release\PoissonRecon.exe --in \data\PointSets\eagle_cleaned.ply --color 16 --depth 10 --out eagle.d.ply --density --voxel eagel.d.iso --voxelDepth 8 --dirichlet
\Research\PoissonRecon\PoissonRecon\Bin\x64\Release\PoissonRecon.exe --in \data\PointSets\eagle_cleaned.ply --color 16 --depth 10 --out eagle.n.ply --density --voxel eagel.n.iso --voxelDepth 8
\Research\PoissonRecon\PoissonRecon\Bin\x64\Release\SurfaceTrimmer.exe --in eagle.d.ply --out eagle.d.trim.ply --trim 6
\Research\PoissonRecon\PoissonRecon\Bin\x64\Release\SurfaceTrimmer.exe --in eagle.n.ply --out eagle.n.trim.ply --trim 6

