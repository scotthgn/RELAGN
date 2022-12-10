echo "initpackage relagn lmod_relagn.dat `pwd`/src/fortran_version/relagn_dir" | xspec
echo "initpackage relqso lmod_relqso.dat `pwd`/src/fortran_version/relqso_dir" | xspec

echo "lmod relagn `pwd`/src/fortran_version/relagn_dir" | xspec
echo "lmod relqso `pwd`/src/fortran_version/relqso_dir" | xspec
