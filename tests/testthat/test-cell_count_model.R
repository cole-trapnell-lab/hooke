context('test-cell_count_model')


test_that('new_cell_count_model works', {

  # On Github Actions
  skip_if_not(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})


test_that('new_cell_count_model works', {

  # Not on Github Actions
  skip_if(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})


test_that('new_cell_count_model problems', {

  set.seed(2016)

})


test_that('select_model works', {

  # On Github Actions
  #  skip_if_not(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  # Not on Github Actions
  skip_if(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})


test_that('select_model works', {

  # On Github Actions
  skip_if_not(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  # Not on Github Actions
  # skip_if(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})


test_that('select_model problems', {

  # On Github Actions
  #  skip_if_not(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  # Not on Github Actions
  #  skip_if(identical(Sys.getenv("GITHUB_ACTIONS"), "true"))

  set.seed(2016)

})

