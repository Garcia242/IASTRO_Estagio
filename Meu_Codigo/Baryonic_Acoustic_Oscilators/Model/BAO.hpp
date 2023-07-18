
// Code generated by stanc v2.31.0
#include <stan/model/model_header.hpp>
namespace BAO_model_namespace {

using stan::model::model_base_crtp;
using namespace stan::math;


stan::math::profile_map profiles__;
static constexpr std::array<const char*, 34> locations_array__ = 
{" (found before start of program)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 39, column 2 to column 10)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 40, column 2 to column 10)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 47, column 4 to column 35)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 49, column 4 to column 26)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 51, column 4 to column 19)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 53, column 4 to column 27)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 62, column 4 to column 22)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 63, column 4 to column 20)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 64, column 4 to column 21)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 65, column 4 to column 24)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 71, column 4 to column 153)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 67, column 17 to line 73, column 3)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 67, column 2 to line 73, column 3)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 99, column 2 to column 21)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 100, column 2 to column 24)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 103, column 2 to column 30)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 109, column 1 to column 2)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 26, column 2 to column 18)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 27, column 2 to column 19)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 28, column 2 to column 22)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 34, column 2 to column 20)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 35, column 2 to column 19)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 4, column 4 to column 23)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 5, column 4 to column 23)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 8, column 4 to column 39)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 3, column 89 to line 9, column 3)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 12, column 4 to column 23)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 13, column 4 to column 23)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 15, column 4 to column 28)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 16, column 4 to column 22)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 17, column 4 to column 26)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 18, column 4 to column 71)",
 " (in '/Users/guilhermegarcia/Desktop/IASTRO_Estagio/Meu_Codigo/Baryonic_Acoustic_Oscilators/Model/BAO.stan', line 11, column 30 to line 20, column 1)"};

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
struct rs_functor__ {
  template <typename T0__,
            stan::require_all_t<stan::is_stan_scalar<T0__>>* = nullptr>
  stan::promote_args_t<T0__>
  operator()(const std::vector<T0__>& theta, std::ostream* pstream__) const;
};

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
      local_scalar_t__ H0 = DUMMY_VAR__;
      current_statement__ = 23;
      H0 = stan::model::rvalue(theta, "theta", stan::model::index_uni(1));
      local_scalar_t__ Om = DUMMY_VAR__;
      current_statement__ = 24;
      Om = stan::model::rvalue(theta, "theta", stan::model::index_uni(2));
      current_statement__ = 25;
      return (1 /
               stan::math::pow(
                 (((Om * stan::math::pow((1 + x), 3)) + 1) - Om), 0.5));
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
      current_statement__ = 27;
      H0 = stan::model::rvalue(theta, "theta", stan::model::index_uni(1));
      local_scalar_t__ Om = DUMMY_VAR__;
      current_statement__ = 28;
      Om = stan::model::rvalue(theta, "theta", stan::model::index_uni(2));
      local_scalar_t__ wm = DUMMY_VAR__;
      current_statement__ = 29;
      wm = (Om * stan::math::pow((H0 / 100), 2));
      local_scalar_t__ wb = DUMMY_VAR__;
      current_statement__ = 30;
      wb = 0.02226;
      local_scalar_t__ wn = DUMMY_VAR__;
      current_statement__ = 31;
      wn = (0.0107 * 0.06);
      current_statement__ = 32;
      return ((55.154 *
                stan::math::exp((-72.3 * stan::math::pow((wn + 0.0006), 2))))
               /
               (stan::math::pow(wm, 0.25351) * stan::math::pow(wb, 0.12807)));
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
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

template <typename T0__, stan::require_all_t<stan::is_stan_scalar<T0__>>*>
stan::promote_args_t<T0__>
rs_functor__::operator()(const std::vector<T0__>& theta,
                         std::ostream* pstream__)  const
{
  return rs(theta, pstream__);
}

 class BAO_model final : public model_base_crtp<BAO_model> {

 private:
  std::vector<double> z;
  std::vector<double> dv;
  std::vector<double> error;
  std::vector<double> x_r;
  std::vector<int> x_i; 
  
 
 public:
  ~BAO_model() { }
  
  inline std::string model_name() const final { return "BAO_model"; }

  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.31.0", "stancflags = "};
  }
  
  
  BAO_model(stan::io::var_context& context__, unsigned int random_seed__ = 0,
            std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double ;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static constexpr const char* function__ = "BAO_model_namespace::BAO_model";
    (void) function__;  // suppress unused var warning
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    (void) DUMMY_VAR__;  // suppress unused var warning
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 18;
      context__.validate_dims("data initialization","z","double",
           std::vector<size_t>{static_cast<size_t>(6)});
      z = std::vector<double>(6, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 18;
      z = context__.vals_r("z");
      current_statement__ = 19;
      context__.validate_dims("data initialization","dv","double",
           std::vector<size_t>{static_cast<size_t>(6)});
      dv = std::vector<double>(6, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 19;
      dv = context__.vals_r("dv");
      current_statement__ = 20;
      context__.validate_dims("data initialization","error","double",
           std::vector<size_t>{static_cast<size_t>(6)});
      error = 
        std::vector<double>(6, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 20;
      error = context__.vals_r("error");
      current_statement__ = 21;
      x_r = std::vector<double>(0, std::numeric_limits<double>::quiet_NaN());
      
      
      current_statement__ = 22;
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
    static constexpr const char* function__ = "BAO_model_namespace::log_prob";
    (void) function__;  // suppress unused var warning
    
    try {
      local_scalar_t__ H0 = DUMMY_VAR__;
      current_statement__ = 1;
      H0 = in__.template read<local_scalar_t__>();
      local_scalar_t__ Om = DUMMY_VAR__;
      current_statement__ = 2;
      Om = in__.template read<local_scalar_t__>();
      std::vector<local_scalar_t__> theta =
         std::vector<local_scalar_t__>(2, DUMMY_VAR__);
      current_statement__ = 3;
      stan::model::assign(theta, std::vector<local_scalar_t__>{H0, Om},
        "assigning variable theta");
      std::vector<local_scalar_t__> dv_theo =
         std::vector<local_scalar_t__>(6, DUMMY_VAR__);
      local_scalar_t__ rf = DUMMY_VAR__;
      current_statement__ = 5;
      rf = 150;
      local_scalar_t__ c = DUMMY_VAR__;
      current_statement__ = 6;
      c = (2.9979 * stan::math::pow(10, 5));
      std::vector<local_scalar_t__> D_A =
         std::vector<local_scalar_t__>(6, DUMMY_VAR__);
      std::vector<local_scalar_t__> H =
         std::vector<local_scalar_t__>(6, DUMMY_VAR__);
      std::vector<local_scalar_t__> Dv =
         std::vector<local_scalar_t__>(6, DUMMY_VAR__);
      std::vector<local_scalar_t__> ztest =
         std::vector<local_scalar_t__>(6, DUMMY_VAR__);
      current_statement__ = 13;
      for (int i = 1; i <= 6; ++i) {
        current_statement__ = 11;
        stan::model::assign(dv_theo,
          ((rf / rs(theta, pstream__)) *
            stan::math::pow(
              (((stan::model::rvalue(z, "z", stan::model::index_uni(i)) *
                  stan::math::pow((c / H0), 3)) *
                 (1.0 /
                   stan::math::pow(
                     (((Om *
                         stan::math::pow(
                           (1 +
                             stan::model::rvalue(z, "z",
                               stan::model::index_uni(i))), 3)) + 1) - Om),
                     (1.0 / 2)))) *
                stan::math::pow(
                  stan::math::integrate_1d(integrand_functor__(), 0,
                    stan::model::rvalue(z, "z", stan::model::index_uni(i)),
                    theta, x_r, x_i, pstream__), 2)), (1.0 / 3))),
          "assigning variable dv_theo", stan::model::index_uni(i));
      }
      {
        current_statement__ = 14;
        lp_accum__.add(stan::math::normal_lpdf<propto__>(H0, 70, 5));
        current_statement__ = 15;
        lp_accum__.add(stan::math::normal_lpdf<propto__>(Om, 0.3, 0.1));
        current_statement__ = 16;
        lp_accum__.add(stan::math::normal_lpdf<propto__>(dv_theo, dv, error));
        ;
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
    static constexpr const char* function__ = "BAO_model_namespace::write_array";
    (void) function__;  // suppress unused var warning
    
    try {
      double H0 = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 1;
      H0 = in__.template read<local_scalar_t__>();
      double Om = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 2;
      Om = in__.template read<local_scalar_t__>();
      std::vector<double> theta =
         std::vector<double>(2, std::numeric_limits<double>::quiet_NaN());
      std::vector<double> dv_theo =
         std::vector<double>(6, std::numeric_limits<double>::quiet_NaN());
      double rf = std::numeric_limits<double>::quiet_NaN();
      double c = std::numeric_limits<double>::quiet_NaN();
      std::vector<double> D_A =
         std::vector<double>(6, std::numeric_limits<double>::quiet_NaN());
      std::vector<double> H =
         std::vector<double>(6, std::numeric_limits<double>::quiet_NaN());
      std::vector<double> Dv =
         std::vector<double>(6, std::numeric_limits<double>::quiet_NaN());
      std::vector<double> ztest =
         std::vector<double>(6, std::numeric_limits<double>::quiet_NaN());
      out__.write(H0);
      out__.write(Om);
      if (stan::math::logical_negation((stan::math::primitive_value(
            emit_transformed_parameters__) || stan::math::primitive_value(
            emit_generated_quantities__)))) {
        return ;
      } 
      current_statement__ = 3;
      stan::model::assign(theta, std::vector<local_scalar_t__>{H0, Om},
        "assigning variable theta");
      current_statement__ = 5;
      rf = 150;
      current_statement__ = 6;
      c = (2.9979 * stan::math::pow(10, 5));
      current_statement__ = 13;
      for (int i = 1; i <= 6; ++i) {
        current_statement__ = 11;
        stan::model::assign(dv_theo,
          ((rf / rs(theta, pstream__)) *
            stan::math::pow(
              (((stan::model::rvalue(z, "z", stan::model::index_uni(i)) *
                  stan::math::pow((c / H0), 3)) *
                 (1.0 /
                   stan::math::pow(
                     (((Om *
                         stan::math::pow(
                           (1 +
                             stan::model::rvalue(z, "z",
                               stan::model::index_uni(i))), 3)) + 1) - Om),
                     (1.0 / 2)))) *
                stan::math::pow(
                  stan::math::integrate_1d(integrand_functor__(), 0,
                    stan::model::rvalue(z, "z", stan::model::index_uni(i)),
                    theta, x_r, x_i, pstream__), 2)), (1.0 / 3))),
          "assigning variable dv_theo", stan::model::index_uni(i));
      }
      if (emit_transformed_parameters__) {
        out__.write(theta);
        out__.write(dv_theo);
        out__.write(rf);
        out__.write(c);
        out__.write(D_A);
        out__.write(H);
        out__.write(Dv);
        out__.write(ztest);
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
      local_scalar_t__ H0 = DUMMY_VAR__;
      H0 = in__.read<local_scalar_t__>();
      out__.write(H0);
      local_scalar_t__ Om = DUMMY_VAR__;
      Om = in__.read<local_scalar_t__>();
      out__.write(Om);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    } // transform_inits_impl() 
    
  inline void get_param_names(std::vector<std::string>& names__) const {
    
    names__ = std::vector<std::string>{"H0", "Om", "theta", "dv_theo", "rf",
      "c", "D_A", "H", "Dv", "ztest"};
    
    } // get_param_names() 
    
  inline void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{},
      std::vector<size_t>{}, std::vector<size_t>{static_cast<size_t>(2)},
      std::vector<size_t>{static_cast<size_t>(6)}, std::vector<size_t>{
      }, std::vector<size_t>{}, std::vector<size_t>{static_cast<size_t>(6)},
      std::vector<size_t>{static_cast<size_t>(6)},
      std::vector<size_t>{static_cast<size_t>(6)},
      std::vector<size_t>{static_cast<size_t>(6)}};
    
    } // get_dims() 
    
  inline void constrained_param_names(
                                      std::vector<std::string>& param_names__,
                                      bool emit_transformed_parameters__ = true,
                                      bool emit_generated_quantities__ = true) const
    final {
    
    param_names__.emplace_back(std::string() + "H0");
    param_names__.emplace_back(std::string() + "Om");
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "theta" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= 6; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "dv_theo" + '.' + std::to_string(sym1__));
        } 
      }
      param_names__.emplace_back(std::string() + "rf");
      param_names__.emplace_back(std::string() + "c");
      for (int sym1__ = 1; sym1__ <= 6; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "D_A" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= 6; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "H" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= 6; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "Dv" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= 6; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "ztest" + '.' + std::to_string(sym1__));
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
    
    param_names__.emplace_back(std::string() + "H0");
    param_names__.emplace_back(std::string() + "Om");
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= 2; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "theta" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= 6; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "dv_theo" + '.' + std::to_string(sym1__));
        } 
      }
      param_names__.emplace_back(std::string() + "rf");
      param_names__.emplace_back(std::string() + "c");
      for (int sym1__ = 1; sym1__ <= 6; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "D_A" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= 6; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "H" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= 6; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "Dv" + '.' + std::to_string(sym1__));
        } 
      }
      for (int sym1__ = 1; sym1__ <= 6; ++sym1__) {
        {
          param_names__.emplace_back(std::string() + "ztest" + '.' + std::to_string(sym1__));
        } 
      }
    }
    
    if (emit_generated_quantities__) {
      
    }
    
    } // unconstrained_param_names() 
    
  inline std::string get_constrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"H0\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"Om\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"theta\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(2) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"dv_theo\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(6) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"rf\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"c\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"D_A\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(6) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"H\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(6) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"Dv\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(6) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"ztest\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(6) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"}]");
    
    } // get_constrained_sizedtypes() 
    
  inline std::string get_unconstrained_sizedtypes() const {
    
    return std::string("[{\"name\":\"H0\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"Om\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"theta\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(2) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"dv_theo\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(6) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"rf\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"c\",\"type\":{\"name\":\"real\"},\"block\":\"transformed_parameters\"},{\"name\":\"D_A\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(6) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"H\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(6) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"Dv\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(6) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"},{\"name\":\"ztest\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(6) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"transformed_parameters\"}]");
    
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
  (((((((2 + 6) + 1) + 1) + 6) + 6) + 6) + 6);
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
  (((((((2 + 6) + 1) + 1) + 6) + 6) + 6) + 6);
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
     constexpr std::array<const char*, 2> names__{"H0", "Om"};
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
     }; } using stan_model = BAO_model_namespace::BAO_model;

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
  return BAO_model_namespace::profiles__;
}

#endif


