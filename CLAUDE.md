#### Fail-Fast, No Fallbacks
- **No Silent Fallbacks**: Code must fail immediately when expected conditions aren't met. Silent fallback behavior masks bugs and creates unpredictable systems.
- **Explicit Error Messages**: When something goes wrong, stop execution with clear error messages explaining what failed and what was expected.
- **Example**: `raise ValueError(f"Required model {model_name} not found")` instead of falling back to first available model.
