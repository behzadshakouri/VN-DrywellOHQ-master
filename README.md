# VN DrywellOHQ

## About

VN-DrywellOHQ contains C++/Qt tools, OpenHydroQual models, and analysis scripts
for drywell and bioretention research. The repository combines model creation,
spatial field and property generation, electrical-resistivity-tomography (ERT)
comparison, time-series inputs, and post-processing utilities.

The included models and data are research assets. Review site assumptions,
units, boundary conditions, material properties, and calibration before reuse.

## OpenHydroQual resources

| Resource | Link |
|---|---|
| Project website | [openhydroqual.com](https://openhydroqual.com/) |
| Documentation and tutorials | [OpenHydroQual knowledge base](https://openhydroqual.com/knowledge-base/) |
| Searchable examples | [examples.openhydroqual.com](https://examples.openhydroqual.com/) |
| Example source repository | [ArashMassoudieh/OpenHydroQual-Examples](https://github.com/ArashMassoudieh/OpenHydroQual-Examples) |
| Source code and releases | [ArashMassoudieh/OpenHydroQual](https://github.com/ArashMassoudieh/OpenHydroQual) |
| Video tutorials | [OpenHydroQual on YouTube](https://www.youtube.com/@openhydroqual218) |

## Repository contents

- `modelcreator.*`: OpenHydroQual model generation;
- `fieldgenerator.*` and `propertygenerator.*`: spatial field and property
  construction;
- `ERT_comparison.*`: ERT/model comparison utilities;
- `VN_ref.ohq`: reference OpenHydroQual project;
- `VN Bioretention/`: bioretention project materials;
- `Data analysis and plot scripts/`: analysis and documentation workflows; and
- precipitation and ERT conversion files: supporting inputs and preprocessing.

## Building

Use the supplied `Makefile` or open `VN-DrywellOHQ.pro` in a compatible Qt
environment. Configure local paths and OpenHydroQual resources before building.

## Status

This is research software, not a design calculator. Model outputs require
site-specific verification and calibration.
