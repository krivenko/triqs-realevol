time_expr te1("t^2"), te2("t + sin(_pi/2)"), te3("sqrt(9.0) + 1.5");
    
// Check whether the expressions really depend on time
if(te1.is_constant() || te2.is_constant() || !te3.is_constant())
    return EXIT_FAILURE;
    
// Write to a archive
std::stringstream archive_str;
boost::archive::text_oarchive oa(archive_str);
oa << te1; oa << te2; oa << te3;
    
// Read from an archive
boost::archive::text_iarchive ia(archive_str);
time_expr read_expr;
ia >> read_expr; if(read_expr != te1) return EXIT_FAILURE;
ia >> read_expr; if(read_expr != te2) return EXIT_FAILURE;
ia >> read_expr; if(read_expr != te3) return EXIT_FAILURE;
    
// Check correctness of numerical expressions
double T[] = {0, 0.1, 10, 55};
double TE1_res[] = {0, 0.01, 100, 3025};
double TE2_res[] = {1, 1.1, 11, 56};
double TE3_res[] = {4.5, 4.5, 4.5, 4.5};

for(int i = 0; i < 4; ++i){
    if(std::abs(te1(T[i]) - TE1_res[i]) >= 1e-10) return EXIT_FAILURE;
    if(std::abs(te2(T[i]) - TE2_res[i]) >= 1e-10) return EXIT_FAILURE;
    if(std::abs(te3(T[i]) - TE3_res[i]) >= 1e-10) return EXIT_FAILURE;
}

// Unary minus
time_expr mte2 = -te2;
for(int i = 0; i < 4; ++i){
    if(std::abs(mte2(T[i]) - (-TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
}
    
// Addition of expressions
time_expr te1pte2 = te1 + te2;
time_expr te1phalf = te1 + 0.5;
time_expr halfpte2 = 0.5 + te2;

for(int i = 0; i < 4; ++i){
    if(std::abs(te1pte2(T[i]) - (TE1_res[i]+TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
    if(std::abs(te1phalf(T[i]) - (TE1_res[i]+0.5)) >= 1e-10) return EXIT_FAILURE;
    if(std::abs(halfpte2(T[i]) - (0.5+TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
}
    
// Subtraction of expressions
time_expr te1mte2 = te1 - te2;
time_expr te1mhalf = te1 - 0.5;
time_expr halfmte2 = 0.5 - te2;
for(int i = 0; i < 4; ++i){
    if(std::abs(te1mte2(T[i]) - (TE1_res[i]-TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
    if(std::abs(te1mhalf(T[i]) - (TE1_res[i]-0.5)) >= 1e-10) return EXIT_FAILURE;
    if(std::abs(halfmte2(T[i]) - (0.5-TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
}
   
 // Multiplication of expressions
time_expr te1ppte2 = te1 * te2;
time_expr te1pphalf = te1 * 0.5;
time_expr halfppte2 = 0.5 * te2;
for(int i = 0; i < 4; ++i){
    if(std::abs(te1ppte2(T[i]) - (TE1_res[i]*TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
    if(std::abs(te1pphalf(T[i]) - (TE1_res[i]*0.5)) >= 1e-10) return EXIT_FAILURE;
    if(std::abs(halfppte2(T[i]) - (0.5*TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
}
    
// Division of expressions 
time_expr te1dte2 = te1 / te2;
time_expr te1dhalf = te1 / 0.5;
time_expr halfdte2 = 0.5 / te2;
for(int i = 0; i < 4; ++i){
    if(std::abs(te1dte2(T[i]) - (TE1_res[i]/TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
    if(std::abs(te1dhalf(T[i]) - (TE1_res[i]/0.5)) >= 1e-10) return EXIT_FAILURE;
    if(std::abs(halfdte2(T[i]) - (0.5/TE2_res[i])) >= 1e-10) return EXIT_FAILURE;
}
    
// Test try_reduce_to_constant()
time_expr te0t("0*t"), te1t("1*t");
if(te0t.is_constant() || te1t.is_constant()) return EXIT_FAILURE;
uniform_mesh<> m(0,100,101);
try_reduce_to_constant(te0t, m);
try_reduce_to_constant(te1t, m);
if(!te0t.is_constant()) return EXIT_FAILURE;
if(te1t.is_constant()) return EXIT_FAILURE;
