library(readxl)
library(dplyr)
library(tidyverse)

setwd("~/Dropbox (The Francis Crick)/TracerX_Renal_WGS/0_data_info_sheets/")

## Sample list - all samples received CRAMs via AWS transfer

trr_rec_df = read.csv("./L_ids_list_aws_received.txt", header = FALSE)
trr_rec_df = trr_rec_df %>%
  mutate(aws_rec="1")


## Helene reconciliation final data

recon_final_df = read_xlsx("./helene_final_reconciliations/TRR sample reconciliation_20240429.xlsx", 
                           sheet = "Total Rec", .name_repair = "universal",
                           col_types = c("text","text","text","text","text","date","text","text","text"))

recon_final_pass = recon_final_df %>%
  filter(In.latest.RE.Release == "YES")

## Helene recon OL with AWS received data

recon_aws_ol = inner_join(recon_final_pass, trr_rec_df,
                         by = c("LPID"="V1"))


## Helene recon OL with AWS received data OL with Flexylims FB infosheet
trr_fb_df = read.csv("./FB_sheets/flexylims_export_TEST_SEQ_husayn_10_11_2022_01_54.csv", colClasses="character")

trr_fb_df1 = trr_fb_df %>%
  filter(status == "Sequencing Completed")

trr_final_ol_df = inner_join(trr_fb_df1, recon_aws_ol,
                             by = c("sub_sample_id"="Lab.Sample.ID"))

### Total of 877 in v3 Release
### Total of 836 overlap with Flexylims database
### All 836 LPIDs have a corresponding CRAM received on AWS

final_dat_df = trr_final_ol_df %>%
  rowwise() %>% mutate(patient_id = str_split(sample_name, "_")[[1]][1]) 

final_dat_selCols_df = final_dat_df %>%
  select(sample_name_cleaned, patient_id, LPID, Participant.ID, sub_sample_id, sample_id, parent_specimen_id, parent_specimen_name, origin, linked_items, batch, Sample.Type, Delivery.Tracker,aws_rec, In.latest.RE.Release)

write.csv(final_dat_selCols_df, "./trr_wgs_allsamples_infosheet_final_v1.csv", row.names = F)

### 

###

## Some of the LPIDs have more than 1 CRAM file
## list here: /nemo/project/proj-turajlics-wgs/TRACERX_WGS/SOURCE/cram_multi_all_list.txt
## some of it also have Dragen aligned CRAM :P

multi_cram_bs = read.table("/Volumes/proj-turajlics-wgs/TRACERX_WGS/SOURCE/cram_multi_all_list.txt")

colnames(multi_cram_bs) <- c("lpid", "delivery_id", "cram")

multi_cram_mergedInfo = inner_join(final_dat_selCols_df, multi_cram_bs, 
                                  by = c("LPID"="lpid"))

length(unique(multi_cram_mergedInfo$LPID))

# write.csv(multi_cram_mergedInfo, "./trr_wgs_multicram_v1.csv", row.names = F)


### Remove patients with any samples with multiple CRAMs and create a samplesheet

pt_id_toRm = unique(multi_cram_mergedInfo$patient_id)

final_dat_filt = final_dat_selCols_df %>%
  filter(!patient_id %in% pt_id_toRm)

length(unique(final_dat_filt$sample_name_cleaned))

# write.csv(final_dat_filt, "./trr_wgs_multcramRm_infosheet_v1.csv", row.names = F)


### Only patients with any samples with multiple CRAMs and create a samplesheet

final_dat_filt_2 = final_dat_selCols_df %>%
  filter(patient_id %in% pt_id_toRm)

length(unique(final_dat_filt$sample_name_cleaned))

# write.csv(final_dat_filt, "./trr_wgs_multcramRm_infosheet_v1.csv", row.names = F)

### =========================

