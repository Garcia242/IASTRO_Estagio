
// Code generated by stanc v2.31.0
#include <stan/model/model_header.hpp>
namespace snia_model_namespace {

using stan::model::model_base_crtp;
using namespace stan::math;


stan::math::profile_map profiles__;
static constexpr std::array<const char*, 38> locations_array__ = 
{" (found before start of program)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 61, column 2 to column 10)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 63, column 2 to column 14)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 70, column 2 to column 37)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 75, column 2 to column 13)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 76, column 2 to column 13)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 77, column 2 to column 13)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 78, column 2 to column 17)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 82, column 4 to column 103)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 83, column 4 to column 29)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 84, column 4 to column 27)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 85, column 4 to column 21)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 80, column 18 to line 87, column 3)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 80, column 2 to line 87, column 3)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 119, column 2 to column 24)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 121, column 2 to column 25)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 128, column 1 to column 22)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 46, column 2 to column 22)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 47, column 2 to column 20)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 48, column 2 to column 21)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 55, column 2 to column 20)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 56, column 2 to column 19)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 7, column 4 to column 23)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 9, column 4 to column 27)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 12, column 4 to column 92)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 4, column 36 to line 17, column 3)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 22, column 4 to column 23)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 24, column 4 to column 27)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 27, column 4 to column 32)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 20, column 89 to line 28, column 3)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 31, column 4 to column 23)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 32, column 4 to column 23)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 33, column 4 to column 22)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 35, column 4 to column 28)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 36, column 4 to column 22)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 37, column 4 to column 26)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 38, column 4 to column 71)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Model4/Model/snia.stan', line 30, column 32 to line 39, column 1)"};

struct rs_functor__ {
  template <typename T0__,
            stan::require_all_t<stan::is_stan_scalar<T0__>>* = nullptr>
  stan::promote_args_t<T0__>
  operator()(const std::vector<T0__>& theta, std::ostream* pstream__) const;
};
struct E_functor__ {
  template <typename T0__, typename T1__,
            stan::require_all_t<stan::is_stan_scalar<T0__>,
                                stan::is_stan_scalar<T1__>>* = nullptr>
  stan::promote_args_t<T0__, T1__>
  operator()(const T0__& x, const std::vector<T1__>& theta,
             std::ostream* pstream__) const;
};
struct integrand_functor__ {
  template <typename T0__, typename T1__, typename T2__, typename T3__,
            stan::require_all_t<stan::is_stan_scalar<T0__>,
                                stan::is_stan_scalar<T1__>,
                                stan::is_stan_scalar<T2__>,
                                stan::is_stan_scalar<T3__>>* = nullptr>
  stan::promote_args_t<T0__, T1__, T2__, T3__>
  operator()(const T0__& x, const T1__& xc, const std::vector<T2__>& theta,
             const std::vector<T3__>& x_r, const std::vector<int>& x_i,
             std::ostream* pstream__) const;
};

template <typename T0__, typename T1__,
          stan::require_all_t<stan::is_stan_scalar<T0__>,
                              stan::is_stan_scalar<T1__>>* = nullptr>
  stan::promote_args_t<T0__, T1__>
  E(const T0__& x, const std::vector<T1__>& theta, std::ostream* pstream__) {
    using local_scalar_t__ = stan::promote_args_t<T0__, T1__>;
    int current_statement__ = 0; 
    static constexpr bool propto__ = true;
    (void) propto__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    try {
      local_scalar_t__ Om = DUMMY_VAR__;
      current_statement__ = 22;
      Om = stan::model::rvalue(theta, "theta", stan::model::index_uni(1));
      local_scalar_t__ lambda = DUMMY_VAR__;
      current_statement__ = 23;
      lambda = stan::model::rvalue(theta, "theta", stan::model::index_uni(2));
      current_statement__ = 24;
      return stan::math::pow(
               ((((3 / (3 - stan::math::pow(lambda, 2))) * Om) *
                  stan::math::pow((1 + x), 3)) +
                 ((1 - ((3 / (3 - stan::math::pow(lambda, 2))) * Om)) *
                   stan::math::pow((1 + x), stan::math::pow(lambda, 2)))),
               0.5);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    }
template <typename T0__, typename T1__, typename T2__, typename T3__,
          stan::require_all_t<stan::is_stan_scalar<T0__>,
                              stan::is_stan_scalar<T1__>,
                              stan::is_stan_scalar<T2__>,
                              stan::is_stan_scalar<T3__>>* = nullptr>
  stan::promote_args_t<T0__, T1__, T2__, T3__>
  integrand(const T0__& x, const T1__& xc, const std::vector<T2__>& theta,
            const std::vector<T3__>& x_r, const std::vector<int>& x_i,
            std::ostream* pstream__) {
    using local_scalar_t__ = stan::promote_args_t<T0__, T1__, T2__, T3__>;
    int current_statement__ = 0; 
    static constexpr bool propto__ = true;
    (void) propto__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    try {
      local_scalar_t__ Om = DUMMY_VAR__;
      current_statement__ = 26;
      Om = stan::model::rvalue(theta, "theta", stan::model::index_uni(1));
      local_scalar_t__ lambda = DUMMY_VAR__;
      current_statement__ = 27;
      lambda = stan::model::rvalue(theta, "theta", stan::model::index_uni(2));
      current_statement__ = 28;
      return (1 / E(x, std::vector<local_scalar_t__>{Om, lambda}, pstream__));
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    }
template <typename T0__,
          stan::require_all_t<stan::is_stan_scalar<T0__>>* = nullptr>
  stan::promote_args_t<T0__>
  rs(const std::vector<T0__>& theta, std::ostream* pstream__) {
    using local_scalar_t__ = stan::promote_args_t<T0__>;
    int current_statement__ = 0; 
    static constexpr bool propto__ = true;
    (void) propto__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    try {
      local_scalar_t__ H0 = DUMMY_VAR__;
      current_statement__ = 30;
      H0 = stan::model::rvalue(theta, "theta", stan::model::index_uni(1));
      local_scalar_t__ Om = DUMMY_VAR__;
      current_statement__ = 31;
      Om = stan::model::rvalue(theta, "theta", stan::model::index_uni(2));
      local_scalar_t__ M = DUMMY_VAR__;
      current_statement__ = 32;
      M = stan::model::rvalue(theta, "theta", stan::model::index_uni(3));
      local_scalar_t__ wm = DUMMY_VAR__;
      current_statement__ = 33;
      wm = (Om * stan::math::pow((H0 / 100), 2));
      local_scalar_t__ wb = DUMMY_VAR__;
      current_statement__ = 34;
      wb = 0.02226;
      local_scalar_t__ wn = DUMMY_VAR__;
      current_statement__ = 35;
      wn = (0.0107 * 0.06);
      current_statement__ = 36;
      return ((55.154 *
                stan::math::exp((-72.3 * stan::math::pow((wn + 0.0006), 2))))
               /
               (stan::math::pow(wm, 0.25351) * stan::math::pow(wb, 0.12807)));
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    }
template <typename T0__, stan::require_all_t<stan::is_stan_scalar<T0__>>*>
stan::promote_args_t<T0__>
rs_functor__::operator()(const std::vector<T0__>& theta,
                         std::ostream* pstream__)  const
{
  return rs(theta, pstream__);
}

template <typename T0__, typename T1__,
          stan::require_all_t<stan::is_stan_scalar<T0__>,
                              stan::is_stan_scalar<T1__>>*>
stan::promote_args_t<T0__, T1__>
E_functor__::operator()(const T0__& x, const std::vector<T1__>& theta,
                        std::ostream* pstream__)  const
{
  return E(x, theta, pstream__);
}

template <typename T0__, typename T1__, typename T2__, typename T3__,
          stan::require_all_t<stan::is_stan_scalar<T0__>,
                              stan::is_stan_scalar<T1__>,
                              stan::is_stan_scalar<T2__>,
                              stan::is_stan_scalar<T3__>>*>
stan::promote_args_t<T0__, T1__, T2__, T3__>
integrand_functor__::operator()(const T0__& x, const T1__& xc,
                                const std::vector<T2__>& theta,
                                const std::vector<T3__>& x_r,
                                const std::vector<int>& x_i,
                                std::ostream* pstream__)  const
{
  return integrand(x, xc, theta, x_r, x_i, pstream__);
}

 class snia_model final : public model_base_crtp<snia_model> {

 private:
  std::vector<double> zcmb;
  std::vector<double> mb;
  std::vector<double> dmb;
  std::vector<double> x_r;
  std::vector<int> x_i; 
  
 
 public:
  ~snia_model() { }
  
  inline std::string model_name() const final { return "snia_model"; }

  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.31.0", "stancflags = "};
  }
  
  
  snia_model(stan::io::var_context& context__,
             unsigned int random_seed__ = 0,
             std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static constexpr const char* function__ = "snia_model_namespace::snia_model";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 17;
      context__.validate_dims("data initialization","zcmb","double",
           std::vector<size_t>{static_cast<size_t>(40)});
      zcmb = 
        std::vector<double>(40, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 17;
      zcmb = context__.vals_r("zcmb");
      current_statement__ = 18;
      context__.validate_dims("data initialization","mb","double",
           std::vector<size_t>{static_cast<size_t>(40)});
      mb = std::vector<double>(40, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 18;
      mb = context__.vals_r("mb");
      current_statement__ = 19;
      context__.validate_dims("data initialization","dmb","double",
           std::vector<size_t>{static_cast<size_t>(40)});
      dmb = 
        std::vector<double>(40, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 19;
      dmb = context__.vals_r("dmb");
      current_statement__ = 20;
      x_r = std::vector<double>(0, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 21;
      x_i = std::vector<int>(0, std::numeric_limits<int>::min());
      
      
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = 1 + 1;
    
  }
  
  template <bool propto__, bool jacobian__ , typename VecR, typename VecI, 
  stan::require_vector_like_t<VecR>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr> 
  inline stan::scalar_type_t<VecR> log_prob_impl(VecR& params_r__,
                                                 VecI& params_i__,
                                                 std::ostream* pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    static constexpr const char* function__ = "snia_model_namespace::log_prob";
    (void) function__;  // suppress unused var warning
    
    try {
      local_scalar_t__ Om = DUMMY_VAR__;
      current_statement__ = 1;
      Om = in__.template read<local_scalar_t__>();
      local_scalar_t__ lambda = DUMMY_VAR__;
      current_statement__ = 2;
      lambda = in__.template read<local_scalar_t__>();
      std::vector<local_scalar_t__> theta =
         std::vector<local_scalar_t__>(2, DUMMY_VAR__);
      current_statement__ = 3;
      stan::model::assign(theta, std::vector<local_scalar_t__>{Om, lambda},
        "assigning variable theta");
      local_scalar_t__ A = DUMMY_VAR__;
      current_statement__ = 4;
      A = 0;
      local_scalar_t__ B = DUMMY_VAR__;
      current_statement__ = 5;
      B = 0;
      local_scalar_t__ C = DUMMY_VAR__;
      current_statement__ = 6;
      C = 0;
      std::vector<local_scalar_t__> Delta =
         std::vector<local_scalar_t__>(40, DUMMY_VAR__);
      current_statement__ = 13;
      for (int i = 1; i <= 40; ++i) {
        current_statement__ = 8;
        stan::model::assign(Delta,
          (stan::model::rvalue(mb, "mb", stan::model::index_uni(i)) -
            (5.0 *
              stan::math::log10(
                ((1.0 +
                   stan::model::rvalue(zcmb, "zcmb",
                     stan::model::index_uni(i))) *
                  stan::math::integrate_1d(integrand_functor__(), 0,
                    stan::model::rvalue(zcmb, "zcmb",
                      stan::model::index_uni(i)), theta, x_r, x_i, pstream__))))),
          "assigning variable Delta", stan::model::index_uni(i));
        current_statement__ = 9;
        A = (A +
              stan::math::pow(
                (stan::model::rvalue(Delta, "Delta",
                   stan::model::index_uni(i)) /
                  stan::model::rvalue(dmb, "dmb", stan::model::index_uni(i))),
                2));
        current_statement__ = 10;
        B = (B +
              (stan::model::rvalue(Delta, "Delta", stan::model::index_uni(i))
                /
                stan::math::pow(
                  stan::model::rvalue(dmb, "dmb", stan::model::index_uni(i)),
                  2)));
        current_statement__ = 11;
        C = (C +
              stan::math::pow(
                stan::model::rvalue(dmb, "dmb", stan::model::index_uni(i)),
                -2));
      }
      {
        current_statement__ = 14;
        lp_accum__.add(stan::math::normal_lpdf<propto__>(Om, 0.3, 0.1));
        current_statement__ = 15;
        lp_accum__.add(stan::math::normal_lpdf<propto__>(lambda, 0, 10));
        current_statement__ = 16;
        lp_accum__.add((-A + (stan::math::pow(B, 2) / C)));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob_impl() 
    
  template <typename RNG, typename VecR, typename VecI, typename VecVar, 
  stan::require_vector_like_vt<std::is_floating_point, VecR>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr, 
  stan::require_vector_vt<std::is_floating_point, VecVar>* = nullptr> 
  inline void write_array_impl(RNG& base_rng__, VecR& params_r__,
                               VecI& params_i__, VecVar& vars__,
                               const bool emit_transformed_parameters__ = true,
                               const bool emit_generated_quantities__ = true,
                               std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    (void) propto__;
    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    int current_statement__ = 0; 
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    constexpr bool jacobian__ = false;
    (void) DUMMY_VAR__;  // suppress unused var warning
    static constexpr const char* function__ = "snia_model_namespace::write_array";
    (void) function__;  // suppress unused var warning
    
    try {
      double Om = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 1;
      Om = in__.template read<local_scalar_t__>();
      double lambda = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 2;
      lambda = in__.template read<local_scalar_t__>();
      std::vector<double> theta =
         std::vector<double>(2, std::numeric_limits<double>::quiet_NaN());
      double A = std::numeric_limits<double>::quiet_NaN();
      double B = std::numeric_limits<double>::quiet_NaN();
      double C = std::numeric_limits<double>::quiet_NaN();
      std::vector<double> Delta =
         std::vector<double>(40, std::numeric_limits<double>::quiet_NaN());
      out__.write(Om);
      out__.write(lambda);
      if (stan::math::logical_negation((stan::math::primitive_value(
            emit_transformed_parameters__) || stan::math::primitive_value(
            emit_generated_quantities__)))) {
        return ;
      } 
      current_statement__ = 3;
      stan::model::assign(theta, std::vector<local_scalar_t__>{Om, lambda},
        "assigning variable theta");
      current_statement__ = 4;
      A = 0;
      current_statement__ = 5;
      B = 0;
      current_statement__ = 6;
      C = 0;
      current_statement__ = 13;
      for (int i = 1; i <= 40; ++i) {
        current_statement__ = 8;
        stan::model::assign(Delta,
          (stan::model::rvalue(mb, "mb", stan::model::index_uni(i)) -
            (5.0 *
              stan::math::log10(
                ((1.0 +
                   stan::model::rvalue(zcmb, "zcmb",
                     stan::model::index_uni(i))) *
                  stan::math::integrate_1d(integrand_functor__(), 0,
                    stan::model::rvalue(zcmb, "zcmb",
                      stan::model::index_uni(i)), theta, x_r, x_i, pstream__))))),
          "assigning variable Delta", stan::model::index_uni(i));
        current_statement__ = 9;
        A = (A +
              stan::math::pow(
                (stan::model::rvalue(Delta, "Delta",
                   stan::model::index_uni(i)) /
                  stan::model::rvalue(dmb, "dmb", stan::model::index_uni(i))),
                2));
        current_statement__ = 10;
        B = (B +
              (stan::model::rvalue(Delta, "Delta", stan::model::index_uni(i))
                /
                stan::math::pow(
                  stan::model::rvalue(dmb, "dmb", stan::model::index_uni(i)),
                  2)));
        current_statement__ = 11;
        C = (C +
              stan::math::pow(
                stan::model::rvalue(dmb, "dmb", stan::model::index_uni(i)),
                -2));
      }
      if (emit_transformed_parameters__) {
        out__.write(theta);
        out__.write(A);
        out__.write(B);
        out__.write(C);
        out__.write(Delta);
      } 
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      } 
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    } // write_array_impl() 
    
  template <typename VecVar, typename VecI, 
  stan::require_vector_t<VecVar>* = nullptr, 
  stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr> 
  inline void transform_inits_impl(VecVar& params_r__, VecI& params_i__,
                                   VecVar& vars__,
                                   std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      local_scalar_t__ Om = DUMMY_VAR__;
      Om = in__.read<local_scalar_t__>();
      out__.write(Om);
      local_scalar_t__ lambda = DUMMY_VAR__;
      lambda = in__.read<local_scalar_t__>();
      out__.write(lambda);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__ = std::vector<std::string>{"Om", "lambda", "theta", "A", "B",
      "C", "Delta"};
    
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{},
      std::vector<size_t>{}, std::vector<size_t>{static_cast<size_t>(2)},
      std::vector<size_t>{}, std::vector<size_t>{}, std::vector<size_t>{
      }, std::vector<size_t>{static_cast<size_t>(40)}};
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "Om");
    param_names__.emplace_back(std::string() + "lambda");
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "theta" + '.' + std::to_string(sym1__));
        } 
      }
      param_names__.emplace_back(std::string() + "A");
      param_names__.emplace_back(std::string() + "B");
      param_names__.emplace_back(std::string() + "C");
      for (int sym1__ = 1; sym1__ <= 40; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "Delta" + '.' + std::to_string(sym1__));
        } 
      }
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // constrained_param_names() 
    
  inline void unconstrained_param_names(
                                        std::vector<std::string>& param_names__,
                                        bool emit_transformed_parameters__ = true,
                                        bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "Om");
    param_names__.emplace_back(std::string() + "lambda");
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "theta" + '.' + std::to_string(sym1__));
        } 
      }
      param_names__.emplace_back(std::string() + "A");
      param_names__.emplace_back(std::string() + "B");
      param_names__.emplace_back(std::string() + "C");
      for (int sym1__ = 1; sym1__ <= 40; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "Delta" + '.' + std::to_string(sym1__));
        } 
      }
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"Om\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"lambda\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"theta\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(2) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"A\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"B\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"C\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"Delta\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(40) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"}]");
    
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"Om\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"lambda\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"theta\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(2) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"A\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"B\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"C\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"Delta\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(40) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"}]");
    
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    inline void write_array(RNG& base_rng,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                            Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                            const bool emit_transformed_parameters = true,
                            const bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      const size_t num_params__ = 
  (1 + 1);
      const size_t num_transformed = emit_transformed_parameters * 
  ((((2 + 1) + 1) + 1) + 40);
      const size_t num_gen_quantities = emit_generated_quantities * 0;
      const size_t num_to_write = num_params__ + num_transformed +
        num_gen_quantities;
      std::vector<int> params_i;
      vars = Eigen::Matrix<double, Eigen::Dynamic, 1>::Constant(num_to_write,
        std::numeric_limits<double>::quiet_NaN());
      write_array_impl(base_rng, params_r, params_i, vars,
        emit_transformed_parameters, emit_generated_quantities, pstream);
    }

    template <typename RNG>
    inline void write_array(RNG& base_rng, std::vector<double>& params_r,
                            std::vector<int>& params_i,
                            std::vector<double>& vars,
                            bool emit_transformed_parameters = true,
                            bool emit_generated_quantities = true,
                            std::ostream* pstream = nullptr) const {
      const size_t num_params__ = 
  (1 + 1);
      const size_t num_transformed = emit_transformed_parameters * 
  ((((2 + 1) + 1) + 1) + 40);
      const size_t num_gen_quantities = emit_generated_quantities * 0;
      const size_t num_to_write = num_params__ + num_transformed +
        num_gen_quantities;
      vars = std::vector<double>(num_to_write,
        std::numeric_limits<double>::quiet_NaN());
      write_array_impl(base_rng, params_r, params_i, vars,
        emit_transformed_parameters, emit_generated_quantities, pstream);
    }

    template <bool propto__, bool jacobian__, typename T_>
    inline T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
                       std::ostream* pstream = nullptr) const {
      Eigen::Matrix<int, -1, 1> params_i;
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }

    template <bool propto__, bool jacobian__, typename T__>
    inline T__ log_prob(std::vector<T__>& params_r,
                        std::vector<int>& params_i,
                        std::ostream* pstream = nullptr) const {
      return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
    }


    inline void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream = nullptr) const final {
      std::vector<double> params_r_vec(params_r.size());
      std::vector<int> params_i;
      transform_inits(context, params_i, params_r_vec, pstream);
      params_r = Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1>>(
        params_r_vec.data(), params_r_vec.size());
    }

  inline void transform_inits(const stan::io::var_context& context,
                              std::vector<int>& params_i,
                              std::vector<double>& vars,
                              std::ostream* pstream__ = nullptr) const {
     constexpr std::array<const char*, 2> names__{"Om", "lambda"};
      const std::array<Eigen::Index, 2> constrain_param_sizes__{1, 1};
      const auto num_constrained_params__ = std::accumulate(
        constrain_param_sizes__.begin(), constrain_param_sizes__.end(), 0);
    
     std::vector<double> params_r_flat__(num_constrained_params__);
     Eigen::Index size_iter__ = 0;
     Eigen::Index flat_iter__ = 0;
     for (auto&& param_name__ : names__) {
       const auto param_vec__ = context.vals_r(param_name__);
       for (Eigen::Index i = 0; i < constrain_param_sizes__[size_iter__]; ++i) {
         params_r_flat__[flat_iter__] = param_vec__[i];
         ++flat_iter__;
       }
       ++size_iter__;
     }
     vars.resize(num_params_r__);
     transform_inits_impl(params_r_flat__, params_i, vars, pstream__);
    } // transform_inits() 
     }; } using stan_model = snia_model_namespace::snia_model;

#ifndef USING_R

// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}

stan::math::profile_map& get_stan_profile_data() {
  return snia_model_namespace::profiles__;
}

#endif


