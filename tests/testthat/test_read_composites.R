# Tests for read_composites

test_that('parse_composite_info works', {

  # If info is not XML, return an empty list
  expect_equal(parse_composite_info(''), list())
  expect_equal(parse_composite_info('CD8 CD68 CK'), list())

  info = '<?xml version="1.0" encoding="utf-16"?>
<ImageDescription>
 <Version>1</Version>
 <PhenotypeSchema>
  <Name>Unnamed</Name>
  <Phenotypes>
   <Entry>
    <Name>CD8+</Name>
    <Color>255,255,255,0</Color>
   </Entry>
   <Entry>
    <Name>CD68+</Name>
    <Color>255,255,0,0</Color>
   </Entry>
   <Entry>
    <Name>other</Name>
    <Color>255,0,0,255</Color>
   </Entry>
  </Phenotypes>
 </PhenotypeSchema>
 <Composite>
  <Name>CD8 CD68 CK</Name>
  <Components>
   <Entry>
    <Marker>DAPI</Marker>
    <Color>255,0,0,255</Color>
   </Entry>
   <Entry>
    <Marker>CD8</Marker>
    <Color>255,255,255,0</Color>
   </Entry>
  </Components>
 </Composite>
</ImageDescription>'

  value = parse_composite_info(info)
  expect_equal(length(value), 4)
  expect_equal(value$composite_name, "CD8 CD68 CK")
  expect_equal(value$components, c(DAPI = "#0000FF", CD8 = "#FFFF00"))
  expect_equal(value$scheme_name, "Unnamed")
  expect_equal(value$phenotypes, c(`CD8+` = "#FFFF00",
                                    `CD68+` = "#FF0000", other = "#0000FF"))
})


test_that('read_composite works', {
  multi_composite_path =
    "C:/Research/phenoptrTestData/Multi-schema/Multi_composite_image.tif"
  skip_if_not(file.exists(multi_composite_path))

  imgs = read_composites(multi_composite_path)
  expect_equal(length(imgs), 3)
  expect_true(all(purrr::map_int(imgs, purrr::attr_getter('width'))==1860))
  expect_true(all(purrr::map_int(imgs, purrr::attr_getter('length'))==1396))

  expected_composite_names = c("Standard", "CD8 FoxP3", "PD-1")
  expect_equal(purrr::map_chr(imgs, purrr::attr_getter('composite_name')),
               expected_composite_names)

  component_names = purrr::map(imgs, purrr::attr_getter('components')) %>%
    purrr::map(names)
  expected_components = list(
    c("CD8", "PD-L1", "FoxP3", "PD-1", "CK", "DAPI", "CD68", "Autofluorescence"),
    c("CD8", "FoxP3", "CK", "DAPI"),
    c("PD-1", "CK", "DAPI"))
  expect_equal(component_names, expected_components)

  # Try reading just the info
  infos = readTIFFDirectory(multi_composite_path, all=TRUE)

  # infos should be a list of lists
  expect_equal(length(infos), 3)
  expect_true(all(purrr::map_lgl(infos, ~inherits(.x, 'list'))))

  # Image tags are now list elements, not attributes
  expect_true(all(purrr::map_int(infos, 'width')==1860))
  expect_true(all(purrr::map_int(infos, 'length')==1396))

  # The description has to be parsed
  parsed_description = infos %>%
    purrr::map('description') %>%
    purrr::map(parse_composite_info)
  expect_equal(purrr::map_chr(parsed_description, 'composite_name'),
               expected_composite_names)

  component_names = purrr::map(parsed_description, 'components') %>%
    purrr::map(names)
  expect_equal(component_names, expected_components)
})

