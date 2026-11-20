# adamantine-utils
Python package with Adamantine-related utilities

## Utilities

- `melt_pool_analysis`: Extracts melt-pool dimensions from Adamantine VTK output and computes summary statistics or plots their time series.
- `validate_input_database`: Validates an Adamantine JSON input database against the expected physics, boundary, and discretization settings.
- `martensite_score`: Scores an Adamantine run by the time each point spends in a material's precipitation temperature band after forming martensite. Includes a 17-4PH preset; other materials are set with command-line temperatures.
