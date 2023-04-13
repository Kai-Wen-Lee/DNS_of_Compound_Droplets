void output_vtu_bin_foreach (scalar * list, vector * vlist, FILE * fp, bool linear)
{
#if defined(_OPENMP)
  int num_omp = omp_get_max_threads();
  omp_set_num_threads(1);
#endif

  vertex scalar marker[];
  int no_points = 0, no_cells=0 ;
  foreach_vertex(){
    marker[] = _k;
    no_points += 1;
  }
  foreach(){
    no_cells += 1;
  }

  fputs ("<?xml version=\"1.0\"?>\n"
  "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n", fp);
  fputs ("\t <UnstructuredGrid>\n", fp);
  fprintf (fp,"\t\t <Piece NumberOfPoints=\"%d\" NumberOfCells=\"%d\">\n", no_points, no_cells);
  fputs ("\t\t\t <CellData Scalars=\"scalars\">\n", fp);

  int count = 0;
  for (scalar s in list) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"%s\" format=\"appended\" offset=\"%d\">\n", s.name,count);
    count += ((no_cells)+1)*8;
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  for (vector v in vlist) {
    fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" Name=\"Vect-%s\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n", v.x.name,count);
    count += ((no_cells*3)+1)*8;
    fputs ("\t\t\t\t </DataArray>\n", fp);
  }
  fputs ("\t\t\t </CellData>\n", fp);
  fputs ("\t\t\t <Points>\n", fp);
  fprintf (fp,"\t\t\t\t <DataArray type=\"Float64\" NumberOfComponents=\"3\"  format=\"appended\" offset=\"%d\">\n",count);
  count += ((no_points*3)+1)*8;
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Points>\n", fp);
  fputs ("\t\t\t <Cells>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n", fp);
  foreach(){
#if dimension == 2
// Edit OLAND:
    // %g will turn into scientific notation for high integer values. this does not work with paraview
    //fprintf (fp, "%g %g %g %g \n", marker[], marker[1,0], marker[1,1], marker[0,1]);
    int ape1 = marker[];
    int ape2 = marker[1,0];
    int ape3 = marker[1,1];
    int ape4 = marker[0,1];
    fprintf (fp, "%u %u %u %u \n", ape1, ape2, ape3, ape4);
#endif
#if dimension == 3
    // Edit OLAND:
    // %g will turn into scientific notation for high integer values. this does not work with paraview
    //fprintf (fp, "%g %g %g %g %g %g %g %g\n", marker[], marker[1,0,0], marker[1,1,0], marker[0,1,0],marker[0,0,1], marker[1,0,1], marker[1,1,1], marker[0,1,1]);
    int ape1 = marker[];
    int ape2 = marker[1,0,0];
    int ape3 = marker[1,1,0];
    int ape4 = marker[0,1,0];
    int ape5 = marker[0,0,1];
    int ape6 = marker[1,0,1];
    int ape7 = marker[1,1,1];
    int ape8 = marker[0,1,1];
    fprintf (fp, "%u %u %u %u %u %u %u %u\n", ape1, ape2, ape3, ape4, ape5, ape6, ape7, ape8);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n", fp);
  for (int i = 1; i < no_cells+1; i++){
#if dimension == 2
    fprintf (fp, "%d \n", i*4);
#endif
#if dimension == 3
    fprintf (fp, "%d \n", i*8);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t\t <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n", fp);
  foreach(){
#if dimension == 2
    fputs ("9 \n", fp);
#endif
#if dimension == 3
    fputs ("12 \n", fp);
#endif
  }
  fputs ("\t\t\t\t </DataArray>\n", fp);
  fputs ("\t\t\t </Cells>\n", fp);
  fputs ("\t\t </Piece>\n", fp);
  fputs ("\t </UnstructuredGrid>\n", fp);
  fputs ("\t <AppendedData encoding=\"raw\">\n", fp);
  fputs ("_", fp);
  unsigned long long block_len=no_cells*8;
#if dimension == 2
  double vz=0;
#endif
  for (scalar s in list) {
    fwrite (&block_len, sizeof (unsigned long long), 1, fp);
    foreach()
      fwrite (&val(s), sizeof (double), 1, fp);
  }
  block_len=no_cells*8*3;
  for (vector v in vlist) {
    fwrite (&block_len, sizeof (unsigned long long), 1, fp);
    foreach(){
      fwrite (&val(v.x), sizeof (double), 1, fp);
      fwrite (&val(v.y), sizeof (double), 1, fp);
#if dimension == 2
      fwrite (&vz, sizeof (double), 1, fp);
#endif
#if dimension == 3
      fwrite (&val(v.z), sizeof (double), 1, fp);
#endif
    }
  }
  block_len=no_points*8*3;
  fwrite (&block_len, sizeof (unsigned long long), 1, fp);
  foreach_vertex(){
    fwrite (&x, sizeof (double), 1, fp);
    fwrite (&y, sizeof (double), 1, fp);
    fwrite (&z, sizeof (double), 1, fp);
  }
  fputs ("\t\n", fp);
  fputs ("\t </AppendedData>\n", fp);
  fputs ("</VTKFile>\n", fp);
  fflush (fp);
#if defined(_OPENMP)
  omp_set_num_threads(num_omp);
#endif
}