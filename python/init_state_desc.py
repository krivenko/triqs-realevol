from wrap_generator import *

# The module
module = module_(full_name = "init_state", app_name = "realevol", doc = "Functions to produce initial states for real-time evolution")

module.use_module('operators', 'realevol')

module.add_include("init_state.hpp")
module.add_include("<triqs/python_tools/converters/set.hpp>")
module.add_include("<triqs/python_tools/converters/map.hpp>")
module.add_using("namespace realevol")
module.add_using("namespace realevol::hilbert_space")

# The class solver
c = class_(
        py_type = "InitState",  # name of the python class
        c_type = "init_state",   # name of the C++ class
        is_printable = True,
        hdf5 = True,
        doc = "Initial state, including information about the Hilbert space structure"
)

module.add_class(c)

module.add_function(name = "make_pure_init_state",
                    signature = """init_state(operator_t generator, std::set<indices_t> fermion_indices, std::set<indices_t> boson_indices = {},                                   std::map<operators::indices_t, int> bits_per_boson = {})""",
                    calling_pattern = "auto result = make_pure_init_state(generator, fundamental_operator_set(fermion_indices,boson_indices), bits_per_boson);")

# TODO make_zerotemp_init_state
# TODO make_thermal_init_state

module.generate_code()
