# This script countains functions to generate HTML tables for variant pairs
library(glue)

# Test data
# vp=read_tsv("variant_lookups/OGI581_vp_phase.tsv")

POP_NAMES = c( 
  'all' = "All samples",
  'afr' = "African/African American",
  'ami' = 'Amish',
  'amr' = "Admixed American",
  'asj' = "Ashkenazi Jewish",
  'eas' = "East Asian",
  'fin' = "Finnish",
  'nfe' = "Non-Finnish European",
  'sas' = "South Asian",
  'oth' = "Other (population not assigned)"
)

read_vp_tsv<-function(path){
  return(
    read_tsv(
      path,
      col_types = cols(
        chrom = col_character(),
        ref1 = col_character(),
        alt1 = col_character(),
        ref2 = col_character(),
        alt2 = col_character()
      )
    )
  )
}

get_vp_markdown <- function(data){
  return(
    data %>%
      add_vp_html_table() %>%
      mutate(pop=factor(
        pop, levels=c(names(POP_NAMES))
        ),
        pop_name = POP_NAMES[pop]
        ) %>%
      arrange(
        chrom,
        pos1,
        pos2,
        pop
      ) %>%
      glue_data(
        "
    ### {pop_name}
    
    
    {html_table}
    
    <b>Estimated compound het probability for {pop_name}: {em_p_chet_adj}</b>
    "
      ) %>%
      paste(collapse = "\n\n\n")
  )
}

add_vp_html_table<-function(data){
  return(
    data %>%
      mutate(
        html_table=glue('
        <table>
          <tr>
            <th rowspan=5 style="transform:rotate(270deg)"> {chrom}:{pos1} {ref1}/{alt1}</th>
            <th colspan=4> {chrom}:{pos2} {ref2}/{alt2}</th>
          </tr>
          <tr style="border-bottom: 1px solid black">
            <th style="border-right: 1px solid black">{pop}</th>
            <td>{ref2}/{ref2}</td>
            <td>{ref2}/{alt2}</td>
            <td>{alt2}/{alt2}</td>
          </tr>
          <tr>
            <td style="border-right: 1px solid black">{ref1}/{ref1}</td>
            <td>{adj.ref_ref}</td>
            <td>{adj.ref_het}</td>
            <td>{adj.ref_hom}</td>
          </tr>
          <tr>
            <td style="border-right: 1px solid black">{ref1}/{alt1}</td>
            <td>{adj.het_ref}</td>
            <td>{adj.het_het}</td>
            <td>{adj.het_hom}</td>
          </tr>
          <tr>
            <td style="border-right: 1px solid black">{alt1}/{alt1}</td>
            <td>{adj.hom_ref}</td>
            <td>{adj.hom_het}</td>
            <td>{adj.hom_hom}</td>
          </tr>
        </table>
                        ')
      )
  )
  
}
