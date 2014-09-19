
#
# development makefile 
#
CC = gcc
#CC = cc
#
# bin for executable
#
#BIN = /usr/local/bin
BIN = ~/bin
#BIN = /RAID/home/dave/bin
#
# development flags
#CFLAGS =  -fullwarn 

#R4000 optimisation flags
#CFLAGS = -O2 -mips2 -sopt 

#R10000 optimisation flags
#CFLAGS = -O2 -fullwarn 

#DEC flags
#CFLAGS  = -newc -g -tune host
#CFLAGS  = -O5 -tune host -fast 

# gnu c flags
CFLAGS = -O6 
#
#CFLAGS = -g -mips2 -xansi -wlint -fullwarn -prototypes 
#
COBJS = main.o centre_of_mass.o moments_of_inertia.o move_molecule.o cube_roots.o \
        eigen_vec_3b3.o write_car.o write_atom_data.o put_string.o mat_transform.o \
        unit_vector.o vec_cross.o size_vector.o read_input.o tok_get.o find_kind.o \
        get_int.o get_integer.o setup_defaults.o read_field.o read_line.o \
        locate_string.o read_hist_info.o int_to_string.o read_hist_frame.o get_doub.o \
        close_contact.o get_double.o atom_separation_squared.o min_image.o \
        cart_to_fract.o fract_to_cart.o cart_latt_vecs.o vec_dot.o mol_dipole.o \
        forster_kappa.o dot_correlation.o read_output.o write_outcsv.o\
        write_csv.o write_distrib_csv.o write_pdb.o generate_neighbours.o standard_bond.o\
        compare_strings.o msd_calc.o read_statis.o gather_molecule.o join_atoms.o \
        radius_gyration.o find_window.o neighbour_order.o circumcircle.o average_rgyr.o std_rgyr.o \
        atomic_bscat_list.o centre_of_bscat.o moments_of_inertia_bscat.o furthest_atom.o \
        build_ellipse.o build_line.o  write_window_csv.o  write_window_csv_titles.o
#
#
LFLAGS = -lm  
#


main.o            : main.c centre_of_mass.c centre_of_bscat.c moments_of_inertia.c eigen_vec_3b3.c write_car.c \
                     dot_correlation.c write_outcsv.c write_pdb.c read_statis.c \
                     mat_transform.c circumcircle.c furthest_atom.c  find_window.c structures.h
read_line.o        : read_line.c maxima.h
locate_string.o    : locate_string.c 
get_doub.o         : get_doub.c get_int.c global_values.h
put_string.o       : put_string.c global_values.h
copy_int.o         : copy_int.c 
next_none_space.o  : next_none_space.c 
int_to_string.o    : int_to_string.c 
next_space.o       : next_space.c 
get_integer.o      : get_integer.c tok_get.c header.h 
get_int.o          : get_int.c maxima.h
generate_neighbours.o : generate_neighbours.c standard_bond.c  atom_separation_squared.c \
                        structures.h maxima.h data.h
standard_bond.o    : standard_bond.c compare_strings.c maxima.h data.h
atom_separation_squared.o : atom_separation_squared.c min_image.c structures.h global_values.h
compare_strings.o  : compare_strings.c 
min_image.o        : min_image.c cart_to_fract.c fract_to_cart.c maxima.h global_values.h structures.h
cart_to_fract.o    : cart_to_fract.c
fract_to_cart.o    : fract_to_cart.c
find_chunk.o       : find_chunk.c maxima.h global_values.h structures.h
move_molecule.o    : move_molecule.c structures.h
rotate_with_flags.o : rotate_with_flags.c move_molecule.c maxima.h global_values.h structures.h
write_car.o        : write_car.c write_atom_data.c structures.h
write_atom_data.o  : write_atom_data.c structures.h
unit_vector.o      : unit_vector.c size_vector.c
size_vector.o      : size_vector.c
move_molecule_with_flags.o : move_molecule_with_flags.c structures.h
string_to_int.o    : string_to_int.c
rotate.o           : rotate.c move_molecule.c global_values.h structures.h
open_file.o        : open_file.c 
print_product.o    : print_product.c
set_unit_matrix.o  : set_unit_matrix.c
make_move.o        : make_move.c unit_vector.c rotate.c size_vector.c \
                     vec_cross.c vec_dot.c rotate_vecs.c \
                     maxima.h global_values.h structures.h constants.h
bra_x_ket.o        : bra_x_ket.c
matrix_x_ket.o      : matrix_x_ket.c
ket_bra_to_matrix.o : ket_bra_to_matrix.c
bra_mat_ket.o       : bra_mat_ket.c
rotate_vecs.o       : rotate_vecs.c maxima.h global_values.h structures.h
abc_from_latt.o     : abc_from_latt.c vec_dot.c vec_cross.c constants.h
get_nonbond_info.o  : get_nonbond_info.c locate_string.c next_none_space.c next_space.c \
                      read_line.c string_from_int.c global_values.h data.h
get_stretch_params.o : get_stretch_params.c locate_string.c string_from_int.c \
                       next_none_space.c next_space.c read_line.c get_doub.c get_int.c \
                       header.h maxima.h structures.h global_values.h data.h
get_equivalences.o   : get_equivalences.c locate_string.c string_from_int.c  \
                       next_none_space.c next_space.c read_line.c get_doub.c get_int.c \
                       header.h maxima.h structures.h global_values.h data.h
string_from_int.o    : string_from_int.c 
setup_defaults.o         : setup_defaults.c maxima.h own_maths.h ewald.h global_values.h \
                           structures.h constants.h data.h header.h
locate_string_in_string.o : locate_string_in_string.c global_values.h 
print_dashes.o            : print_dashes.c header.h 
centre_of_mass.o          : centre_of_mass.c atomic_mass_list.c structures.h
centre_of_bscat.o          : centre_of_bscat.c atomic_bscat_list.c structures.h
atomic_mass_list.o        : atomic_mass_list.c maxima.h data.h
atomic_bscat_list.o        : atomic_bscat_list.c maxima.h data.h
apply_matrix.o            : apply_matrix.c  maxima.h global_values.h structures.h
gen_rot_matrix.o          : gen_rot_matrix.c maxima.h global_values.h structures.h
apply_matrix_to_vecs.o    : apply_matrix_to_vecs.c maxima.h global_values.h structures.h
timer.o                   : timer.c
real_random.o             : real_random.c 
read_input.o              : read_input.c tok_get.c reader.h  maxima.h global_values.h 
find_kind.o               : find_kind.c
etime_me.o                : etime_me.c
calc_det.o                : calc_det.c 
moments_of_inertia.o      : moments_of_inertia.c move_molecule.c cube_roots.c structures.h
moments_of_inertia_bscat.o  : moments_of_inertia_bscat.c move_molecule.c cube_roots.c structures.h
furthest_atom.o           : furthest_atom.c move_molecule.c structures.h
radius_gyration.o         : radius_gyration.c min_image.c structures.h
cube_roots.o              : cube_roots.c own_maths.h 
mat_transform.o           : mat_transform.c
eigen_vec_3b3.o           : eigen_vec_3b3.c unit_vector.c vec_cross.c maxima.h global_values.h \
                            structures.h constants.h
read_input.o              : read_input.c tok_get.c find_kind.c get_double.c 
find_kind.o               : find_kind.c maxima.h reader.h
tok_get.o                 : tok_get.c find_kind.c maxima.h global_values.h reader.h
get_int.o                 : get_int.c global_values.h
read_field.o              : read_field.c locate_string.c find_field.c read_line.c \
                            get_int.c get_doub.c maxima.h global_values.h structures.h \
                            constants.h data.h 
read_hist_info.o          : read_hist_info.c read_line.c int_to_string.c get_int.c \
                            maxima.h global_values.h structures.h constants.h data.h
read_hist_frame.o         : read_hist_frame.c read_line.c int_to_string.c get_int.c \
                            get_doub.c locate_string.c atomic_bscat_list.c maxima.h global_values.h \
                            structures.h constants.h data.h
close_contact.o           : close_contact.c structures.h global_values.h
get_double.o              : get_double.c tok_get.c maxima.h structures.h global_values.h \
                            header.h
cart_latt_vecs.o          : cart_latt_vecs.c vec_dot.c vec_cross.c constants.h own_maths.h \
                            ewald.h
vec_dot.o                 : vec_dot.c 
mol_dipole.o              : mol_dipole.c maxima.h structures.h reader.h header.h
forster_kappa.o           : forster_kappa.c unit_vector.c mol_dipole.c centre_of_mass.c \
                            structures.h
dot_correlation.o         : dot_correlation.c
read_output.o             : read_output.c get_int.c get_doub.c read_line.c locate_string.c \
                            maxima.h global_values.h structures.h constants.h data.h
write_outcsv.o            : write_outcsv.c write_csv.c
write_csv.o               : write_csv.c
write_distrib_csv.o       : write_distrib_csv.c
write_pdb.o               : write_pdb.c maxima.h structures.h
msd_calc.o                : msd_calc.c
read_statis.o             : read_statis.c
gather_molecule.o         : gather_molecule.c
join_atoms.o              : join_atoms.c
find_window.o              : find_window.c
neighbour_order.o          : neighbour_order.c
circumcirlce.o             : circumcircle.c vec_cross.c vec_dot.c unit_vector.c size_vector.c min_image.c \
                             maxima.h global_values.h structures.h
average_rgyr.o             : average_rgyr.c
std_rgyr.o                 : std_rgyr.c
build_ellipse.o            : build_ellipse.c structures.h maxima.h
build_line.o               : build_line.c structures.h maxima.h
write_window_csv.o         :  write_window_csv.c
write_window_csv_titles.o  :  write_window_csv_titles.c

mom_inert : $(COBJS) makefile
		$(CC) $(CFLAGS)  -o mom_inert $(COBJS) $(LFLAGS)

analyse_hist : $(COBJS) makefile
		 $(CC) $(CFLAGS)  -o $(BIN)/analyse_hist $(COBJS) $(LFLAGS)

analyse_hist_test : $(COBJS) makefile
		 $(CC) $(CFLAGS)  -o $(BIN)/analyse_hist_test $(COBJS) $(LFLAGS)

analyse_hist_check : $(COBJS) makefile
		 $(CC) $(CFLAGS)  -o $(BIN)/analyse_hist_check $(COBJS) $(LFLAGS)
