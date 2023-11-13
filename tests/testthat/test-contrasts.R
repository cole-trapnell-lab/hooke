context('test-contrasts')


test_that('estimate_abundances works', {

  # On Github Actions
  skip_if_not(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})


test_that('estimate_abundances works', {

  # Not on Github Actions
  skip_if(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})


test_that('estimate_abundances problems', {

  # On Github Actions
  # skip_if_not(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})


test_that('compare_abundances works', {

  # On Github Actions
  skip_if_not(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})


test_that('compare_abundances works', {

  # Not on Github Actions
  skip_if(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})


test_that('compare_abundances problems', {

  # Not on Github Actions
  # skip_if(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})

