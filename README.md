# Gordon1966 - FEBio Plugin

A C++ FEBio plugin implementing the Gordon et al. (1966) muscle contraction model for finite element biomechanical simulations.

## Overview

**Gordon1966** is a plugin for the FEBio (Finite Element for Biomechanics) framework that extends the standard transversely isotropic Mooney-Rivlin material model with an active contraction stress term based on the force-length relationship described in the seminal Gordon et al. (1966) paper. This plugin is designed for simulating muscle mechanics, specifically capturing the nonlinear behavior of active muscle contraction.

The plugin combines:
- **Mooney-Rivlin hyperelasticity** - For passive material response
- **Fiber-directed anisotropy** - To model directional material properties
- **Gordon force-length curve** - To represent active muscle contraction mechanics

## Features

- **Accurate Muscle Mechanics**: Implements the physiologically validated Gordon et al. (1966) force-length relationship for muscle contraction
- **Transversely Isotropic Material**: Single preferred fiber direction for muscle-aligned stress
- **Hyperelastic Formulation**: Uncoupled (volumetrically incompressible) formulation for better numerical stability
- **FEBio Integration**: Seamless integration with FEBio's powerful finite element framework
- **Cross-Platform**: Pre-compiled binaries for Linux (x64)
- **Extensible Design**: Written in C++ with clean plugin architecture

## Technical Stack

- **Language:** C++
- **Framework:** FEBio (Finite Element for Biomechanics)
- **License:** GPLv3
- **Build System:** Visual Studio project files (included)
- **Platforms:** Linux (x64), compatible with FEBio 3.x SDK

## Repository Structure

```
source/
  FEGordon1966.cpp        Main material class implementation
  FEGordon1966.h          Material class header
  FENewFiberMaterial.cpp  Fiber material with Gordon contraction
  FENewFiberMaterial.h    Fiber material header
  dllmain.cpp             Plugin initialization and registration
  version.h               Version information
  stdafx.h                Precompiled headers
  targetver.h             Target platform version

lib/
  libgordon1966_lnx64.so  Pre-compiled shared library (Linux x64)

example/
  febio.xml               FEBio plugin configuration file
  fix.feb                 Example muscle model file
  fix2.feb                Example muscle model file (variant)
  fix.xplt                Visualization output (binary format)
  fix2.xplt               Visualization output (binary format)
  fix.log                 Simulation log
  fix2.log                Simulation log

Gordon1966.pdf            Reference paper or documentation
LICENSE                   GPLv3 License
README.md                 This file
```

## Material Model

### Mathematical Formulation

The Gordon1966 material model is an uncoupled hyperelastic material that combines:

1. **Isotropic Mooney-Rivlin Part** (passive response):
   - Parameters: `c1`, `c2`
   - Strain energy: W = c1(I₁-3) + c2(I₂-3)

2. **Anisotropic Fiber Part** (passive anisotropy):
   - Parameters: `c3`, `c4`, `c5`, `lam_max`
   - Fiber reinforcement in preferred direction

3. **Active Contraction** (Gordon force-length curve):
   - Activation parameter: `ascl` (activation level, 0-1)
   - Maximum stress: `smax`
   - Force-length curve parameters: `ax`, `ay`, `bx`, `by`, `cx`, `cy`, `dx`, `dy`, `ex`, `ey`

### Material Parameters

| Parameter | Units | Description |
|-----------|-------|-------------|
| `c1` | Pressure | Mooney-Rivlin coefficient for isotropic response |
| `c2` | Pressure | Mooney-Rivlin coefficient for isotropic response |
| `c3` | Pressure | Fiber material coefficient |
| `c4` | Dimensionless | Fiber material exponent |
| `c5` | Pressure | Fiber material coefficient |
| `lam_max` | Dimensionless | Maximum stretch in fiber direction |
| `ascl` | Dimensionless (0-1) | Activation level (muscle contraction magnitude) |
| `smax` | Pressure | Maximum active stress at optimal length |
| `ax`, `ay`, `bx`, `by`, etc. | Dimensionless | Force-length curve fitting parameters |
| `fiber` | Vector | Local fiber direction (typically input as "vector") |

## Installation

### Using Pre-compiled Binary

1. Copy `libgordon1966_lnx64.so` to an accessible location
2. Update your FEBio configuration file (XML) to point to the library:
   ```xml
   <import>/path/to/libgordon1966_lnx64.so</import>
   ```

### Building from Source

#### Prerequisites

- Visual Studio 2019 or later (for Windows)
- C++ compiler with C++11 support (for Linux)
- FEBio SDK 3.x headers and libraries
- CMake or the project build files

#### Build Steps (Linux)

```bash
# With appropriate compiler and FEBio SDK setup:
g++ -shared -fPIC -std=c++11 \
  -I/path/to/febio/sdk/include \
  source/FEGordon1966.cpp \
  source/FENewFiberMaterial.cpp \
  source/dllmain.cpp \
  -L/path/to/febio/sdk/lib -lfebio_core \
  -o libgordon1966_lnx64.so
```

## Usage

### FEBio Configuration

Create or modify a FEBio configuration file (`febio.xml`):

```xml
<?xml version="1.0" encoding="ISO-8859-1"?>
<febio_config version="3.0">
    <default_linear_solver type="pardiso"></default_linear_solver>
    <import>/path/to/libgordon1966_lnx64.so</import>
</febio_config>
```

### FEBio Input File

Define the material in your FEBio model file (`.feb`):

```xml
<material id="1" name="MuscleMaterial" type="gordon1966">
    <c1>1000</c1>                    <!-- Mooney-Rivlin coefficient 1 -->
    <c2>500</c2>                     <!-- Mooney-Rivlin coefficient 2 -->
    <c3>2000</c3>                    <!-- Fiber material coefficient -->
    <c4>2.0</c4>                     <!-- Fiber material exponent -->
    <c5>1500</c5>                    <!-- Fiber material coefficient -->
    <lam_max>1.6</lam_max>           <!-- Maximum stretch -->
    
    <ascl>0.5</ascl>                 <!-- Activation level (0.0 = passive, 1.0 = fully active) -->
    <smax>50000</smax>               <!-- Maximum active stress (Pa) -->
    
    <!-- Gordon force-length curve parameters -->
    <ax>0.8</ax>
    <ay>0.1</ay>
    <bx>1.0</bx>
    <by>0.5</by>
    <cx>1.2</cx>
    <cy>1.0</cy>
    <dx>1.4</dx>
    <dy>0.8</dy>
    <ex>1.6</ex>
    <ey>0.2</ey>
    
    <fiber type="vector">1.0, 0.0, 0.0</fiber>  <!-- Fiber direction -->
</material>
```

### Running a Simulation

```bash
# Run FEBio with the plugin loaded
febio4 -i mymodel.feb -o output.xplt

# Or with explicit config file
FEBIO_CONFIG=febio.xml febio4 -i mymodel.feb
```

## Example Files

The `example/` directory contains sample muscle models:

- **fix.feb** - Example muscle contraction model with specific boundary conditions
- **fix2.feb** - Variant example with different parameters
- **fix.log** / **fix2.log** - Detailed simulation logs with convergence information
- **fix.xplt** / **fix2.xplt** - Binary output files for visualization in PostView or similar tools

These files demonstrate:
- How to set up a muscle tissue model
- Parameter values for realistic muscle behavior
- Boundary conditions and loading scenarios
- Expected output format

## Implementation Details

### Key Classes

#### FEGordon1966
Main material class implementing the composite hyperelastic model:
- Inherits from `FEUncoupledMaterial` for volumetric incompressibility
- Combines isotropic Mooney-Rivlin response with anisotropic fiber contribution
- Methods:
  - `Init()` - Initialize material parameters and check validity
  - `DevStress()` - Calculate deviatoric stress tensor
  - `DevTangent()` - Calculate consistent material stiffness matrix
  - `CreateMaterialPointData()` - Create material point data structures

#### FENewFiberMaterial
Extends fiber material with active contraction:
- Manages fiber-directed anisotropy
- Implements Gordon force-length relationship
- Handles activation-dependent stress generation
- Provides both stress and tangent contributions

### Stress Calculation

The deviatoric stress is calculated as:

```
s = 2/J * dev(T) + s_fiber

where:
  T = B*(W₁ + W₂*I₁) - B₂*W₂
  B = left Cauchy-Green tensor
  W₁, W₂ = strain energy derivatives
  s_fiber = active contraction stress from Gordon curve
```

### Tangent Stiffness

A consistent material tangent is computed for Newton-Raphson iterations:
- Includes fourth-order elasticity tensor contributions
- Accounts for both passive and active material properties
- Ensures numerical stability in implicit FEA

## Force-Length Relationship

The Gordon et al. (1966) force-length curve models the characteristic nonlinear relationship of skeletal muscle:
- **Ascending limb**: Increasing force with stretch (filament overlap increasing)
- **Plateau**: Near-constant maximum force (maximum overlap)
- **Descending limb**: Decreasing force with further stretch (filament overlap decreasing)

The curve parameters (`ax`, `ay`, `bx`, `by`, `cx`, `cy`, `dx`, `dy`, `ex`, `ey`) define a piecewise polynomial approximation of this relationship.

## References

- **Gordon, A. M., Huxley, A. F., & Julian, F. J. (1966).** "The variation in isometric tension with sarcomere length in vertebrate muscle fibres." The Journal of Physiology, 184(1), 170-192.
- **FEBio Documentation:** https://febio.org/
- **FEBio Plugins:** https://febio.org/plugins/

## Troubleshooting

### Plugin Not Loading

**Problem:** FEBio cannot load the plugin library.

**Solutions:**
- Verify the path in `febio.xml` is correct and the file exists
- Check file permissions: `chmod +x libgordon1966_lnx64.so`
- Ensure SDK version compatibility between plugin and FEBio
- Check system architecture matches (x64 vs x32)

### Material Not Recognized

**Problem:** FEBio doesn't recognize the "gordon1966" material type.

**Solutions:**
- Ensure the plugin is loaded before material definition in FEB file
- Check FEBio log for plugin initialization errors
- Verify material type string in XML is exactly "gordon1966"

### Convergence Issues

**Problem:** Simulation diverges with increased error.

**Solutions:**
- Reduce time step size
- Check material parameter ranges (all should be positive)
- Verify activation level is between 0.0 and 1.0
- Ensure fiber direction is properly normalized
- Experiment with different linear solvers (PARDISO, SKYLINE, etc.)

### Physically Unrealistic Results

**Problem:** Stress or displacement values seem incorrect.

**Solutions:**
- Verify material parameter values are appropriate for your tissue type
- Check units consistency (usually SI units: Pa, m)
- Validate the force-length curve parameters against literature values
- Ensure mesh quality and element formulation are suitable

## Contributing

Contributions are welcome! Areas for improvement:
- Additional material models
- Optimization of force-length curve fitting
- Better documentation and examples
- Performance enhancements
- Support for additional platforms (Windows, macOS)

## Performance Considerations

- Plugin performance is largely dependent on FEBio's solver efficiency
- Material point stress/tangent calculations are O(1) complexity
- Pre-compilation to shared library provides good performance
- For large-scale simulations, consider:
  - Mesh coarsening where appropriate
  - Exploiting symmetry
  - Parallel processing if FEBio version supports it

## License

GNU General Public License v3.0 (GPLv3) – see `LICENSE`

This means the code is free and open-source, but derivative works must also be open-source.

## Citation

If you use this plugin in your research, please cite:

```
Stark, H. (2025). Gordon1966 - FEBio Plugin for Muscle Contraction Modeling.
GitHub: https://github.com/heikostark/Gordon1966
```

And reference the original Gordon paper:
```
Gordon, A. M., Huxley, A. F., & Julian, F. J. (1966).
The variation in isometric tension with sarcomere length in vertebrate muscle fibres.
The Journal of Physiology, 184(1), 170-192.
```

## Author

**Heiko Stark** - Research in computational biomechanics and muscle physiology

Contact: [Research Interests](https://stark-jena.de/research-interests/simulation/contraction/material-models/)

## Additional Resources

- **FEBio Official Website:** https://febio.org/
- **FEBio Plugin Registry:** https://febio.org/plugins/
- **Gordon et al. 1966 Paper:** [Available via academic databases]
- **Muscle Physiology References:** Standard biomechanics and physiology textbooks

---

## Development Notes

### Version History

- **v1.0** - Initial release with Gordon force-length curve implementation
- Based on FEBio SDK 3.x
- Compatible with transversely isotropic hyperelastic framework

### Known Limitations

- Currently supports single fiber direction (not multi-directional fibers)
- Activation is uniform across material (no spatial variation)
- No damage or fatigue modeling
- Requires positive-definite material parameters for stability

### Future Enhancements

- Multi-directional fiber families
- Spatially varying activation
- Improved fiber-matrix interaction modeling
- Integration with other biomechanical features
- GPU acceleration for large-scale problems
