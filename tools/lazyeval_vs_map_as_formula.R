# Note to self: This is why `select_rows` uses `as.formula(..., globalenv())`
df = tibble::tibble(x=1:2)

# This fails
f = purrr::map("~x", ~as.formula(.x))
lazyeval::f_eval(f[[1]], df)

# It works if the purrr::map does not use formula notation:
f = purrr::map("~x", as.formula)
lazyeval::f_eval(f[[1]], df)

# Also works if I specify globalenv() in as.formula()
f = purrr::map("~x", ~as.formula(.x, globalenv()))
lazyeval::f_eval(f[[1]], df)

# Or if I make a named function instead of using the formula shortcut
make_formula = function(s) as.formula(s)
f = purrr::map("~x", make_formula)
lazyeval::f_eval(f[[1]], df)

