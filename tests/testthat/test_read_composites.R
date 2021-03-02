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
