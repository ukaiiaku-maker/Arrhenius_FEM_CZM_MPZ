# V2 native-state factory audit

The PF and FEM/CZM factories instantiate their exact source constructors and perform exactly one lossless initialization clone into the V2 lane. All 24 backend/class/temperature fixtures have independent objects, identical complete initial transaction states, and `future_state_imports = 0`. Arrays retain dtype, shape, and values; FEM cleavage/emission RNG streams and thresholds are included.
