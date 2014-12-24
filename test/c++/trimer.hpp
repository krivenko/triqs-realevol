struct write_obs : public boost::static_visitor<> {

    h5::group & gr;
    std::string name;
    write_obs(h5::group & gr, std::string name) : gr(gr), name(name) {}

    template<typename MeshContainer>
    void operator()(MeshContainer const& sol) const { h5_write(gr,name,sol); }
};

void write_results(var_results_t const& res, h5::group gr) {
    for(auto const& r : res)
        boost::apply_visitor(write_obs(gr,r.first), r.second);
}

// Real-valued version of the solver
using solver_t = solver<false>;

int main()
{
    std::array<int,3> atoms {1,2,3};
    std::array<std::string,2> spins {"up","dn"};

    std::set<solver_t::indices_t> all_indices;
    for(auto a : atoms)
        for(auto s : spins)
            all_indices.insert({a,s});

    auto H = operator_t<false>();

    // Chemical potential
    for(auto a : atoms)
        for(auto s : spins)
            H += -mu*n<time_expr_r>(a,s);
    // Hubbard interaction
    for(auto a : atoms) H += U*n<time_expr_r>(a,"up")*n<time_expr_r>(a,"dn");
    // Hopping between atoms
    for(auto s : spins){
        for(auto a1 : atoms)
        for(auto a2 : atoms){
            if(a1==a2) continue;
            H += -V*c_dag<time_expr_r>(a1,s)*c<time_expr_r>(a2,s);
        }
    }

    // Time mesh
    segment_mesh mesh(0,50.0,1001);

    // Solver object
    solver_t S(all_indices);

    // Parameters
    auto params = solve_parameters_t<false>(H,mesh);
    params.verbosity = 2;
    params.stored_psi_values = 20;

    // Observables
    std::map<std::string,operator_t<false>> observables;

    params.observables["unity"] = operator_t<false>() + 1.0; // ugly
    for(auto a : atoms)
        for(auto s : spins) {
            auto name = "n_" + std::to_string(a) + "_" + s;
            params.observables[name] = n<time_expr_r>(a,s);
        }

    params.method = method;

    h5::file output_file(output_filename.c_str(),H5F_ACC_TRUNC);
    h5::group root_gr(output_file);

    // psi0: Sz=1/2
    {
    S.psi0({{1,"up"},{2,"up"},{3,"dn"}}) = 1.0;
    S.solve(params);

    auto res = S.get_results();

    auto gr = root_gr.create_group("Sz_0.5");
    write_results(res,gr);
    }

    // psi0: Sz=0
    {
    S.psi0().amplitudes()() = .0;
    S.psi0({{1,"up"},{1,"dn"},{2,"up"},{3,"dn"}}) = 1.0;
    S.solve(params);

    auto res = S.get_results();

    auto gr = root_gr.create_group("Sz_0");
    write_results(res,gr);
    }

    return EXIT_SUCCESS;
}