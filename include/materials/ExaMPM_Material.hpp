#ifndef EXAMPM_MATERIAL_MATERIAL_H_
#define EXAMPM_MATERIAL_MATERIAL_H_

#include <limits>

#include "Eigen/Dense" // Will eigen be needed at all?
#include "json.hpp"

#include "factory.h"
#include "logger.h"
#include "map.h"
#include "material_utility.h"
#include "particle.h" 
#include "particle_base.h"

// JSON
using Json = nlohmann::json;

namespace ExaMPM {

// Forward declaration of ParticleBase
template <unsigned Tdim>
class ParticleBase;

//! Material base class
//! \brief Base class that stores the information about materials
//! \details Material class stresses and strains
//! \tparam Tdim Dimension

// Data needs to live in ExaMPM_ProblemManager 
template <unsigned Tdim> 
class Material {
 public:
  //! Define a vector of 6 dof
  using Vector6d = Eigen::Matrix<double, 6, 1>;
  //! Define a Matrix of 6 x 6
  using Matrix6x6 = Eigen::Matrix<double, 6, 6>;

  // Constructor with id
  //! \param[in] id Material id
  Material(unsigned id, const Json& material_properties) : id_{id} {
    //! Logger
    std::string logger = "material::" + std::to_string(id);
    console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
  }

  //! Destructor
  virtual ~Material(){};

  //! Delete copy constructor
  Material(const Material&) = delete;

  //! Delete assignement operator
  Material& operator=(const Material&) = delete;

  //! Return id of the material
  unsigned id() const { return id_; }

  //! Get material property
  //! \tparam Ttype Return type for proerpty
  //! \param[in] key Material property key
  //! \retval result Value of material property
  template <typename Ttype>
  Ttype property(const std::string& key);

  //! Initialise history variables
  // Will need to change dense_map function to match what is in ExaMPM
  virtual ExaMPM::dense_map initialise_state_variables() = 0;

  //! State variables
  virtual std::vector<std::string> state_variables() const = 0;

  //! Compute stress
  //! \param[in] stress Stress
  //! \param[in] dstrain Strain
  //! \param[in] particle Constant point to particle base
  //! \param[in] state_vars History-dependent state variables
  //! \retval updated_stress Updated value of stress
  virtual Vector6d compute_stress(const Vector6d& stress,
                                  const Vector6d& dstrain,
                                  const ParticleBase<Tdim>* ptr,
                                  ExaMPM::dense_map* state_vars) = 0;

 
 protected:
  //! material id
  unsigned id_{std::numeric_limits<unsigned>::max()};
  //! Material properties
  Json properties_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // Material class
}  // namespace mpm

#include "ExaMPM_Material.tcc"

#endif  // EXAMPM_MATERIAL_MATERIAL_H_
