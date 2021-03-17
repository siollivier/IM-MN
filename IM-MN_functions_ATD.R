library(data.table)
library(svDialogs)
library(tcltk)
library(OneR)
library(plotly)
library(Rfast)

MZmine_split_list = function(csv){
  
  # We select the intensities by file
  raw=grep('raw filtered Peak height',colnames(csv))
  name_list = list()
  data_list=list()  
  # We initialize the cluster index
  c_index=0
  
  for (i in raw){
    
    # Individualizing file names
    run_name = strsplit(colnames(csv[,i]), split=".raw")[[1]][1]
    
    # Increment cluster index
    c_index=c_index+1
    name_list[c_index]=run_name
    
    # Individualizing each mobilogram from the csv
    a=which(csv[,i]!=0)
    data_list[[c_index]]= subset(csv, subset = csv[,i]>0.002*max(csv[,i]), select=c(1,2,i,i+1)) 
    try_charge = which(data_list[[c_index]][,4]==0)
    data_list[[c_index]][try_charge,4]=1
  }
  
  # merging the name list and data list
  raw_list = list(name_list, data_list) 
  return(raw_list)
}

#######################################################################################
correct_AT_2_DT = function(data_as_list, metadata){
  
  for(i in 1:length((data_as_list[[1]]))){
    
    # looking for the corresponding metadata
    ID = which(metadata$file_name==data_as_list[[1]][[i]])
    
    # extracting dt correction factor
    a=metadata[ID,2] %>% as.numeric()
    data_as_list[[2]][[i]][,2] = data_as_list[[2]][[i]][,2]-a
    
  }
  return(data_as_list)
}

#######################################################################################
CCS_calibration = function(drift_as_list, calib, delay=1.5, gas_mass=28, cat_mass=7.01545){
  
  # creating CCS calibration (log method)
  
  calib$TD = calib$TA-10
  calib$CCSp = (calib$CCSN2/(abs(calib$z)*(1/calib$MW + 1/gas_mass)^0.5))
  calib$DTp = (calib$TD-0.001*delay*calib$mz^0.5)
  calib$ln_CCSp = log(calib$CCSp)
  calib$ln_DTp = log(calib$DTp)
  
  calib_mod=lm(calib$ln_CCSp ~ calib$ln_DTp)
  coef_calib = coef(calib_mod)
  coef_exp = exp(coef_calib[1]) %>% as.numeric()
  slope = coef_calib[2] %>% as.numeric()
  OO_exp = coef_calib[1] %>% as.numeric()
  
  R2 = round(cor(x=calib$ln_DTp, y=calib$ln_CCSp)^2,4)

  calib$ln_CCS_th = slope * calib$ln_DTp + OO_exp
  calib$CCS_mes = calib$DTp^slope*1*(1/calib$MW+1/gas_mass)^0.5
  calib$CCS_mes = round(calib$CCS_mes * coef_exp,2)
  calib$pct_diff_th = round(abs((calib$CCSN2-calib$CCS_mes)/calib$CCSN2*100),2)
  
  mean_dif = round(mean(calib$pct_diff_th),2)
  
  # Plot calibration
  
  ft = list(
    family = "Arial, sans-serif",
    size = 24,
    color = "black"
  )
  
  f1 = list(
    family = "Arial, sans-serif",
    size = 18,
    color = "black"
  )
  
  f2 = list(
    family = "Arial, sans-serif",
    size = 14,
    color = "black"
  )
  
  titre = list(
    text = paste('Log CCS calibration\n R� = ', R2, ' - Mean CCS difference = ', mean_dif,' %',sep=""),
    font = ft
  )
  
  axis_x = list(
    title = "ln (dt')",
    titlefont = f1,
    showticklabels = TRUE,
    tickfont = f2
  )
  
  axis_y = list(
    title = "ln (CCS')",
    titlefont = f1,
    showticklabels = TRUE,
    tickfont = f2
  )
  
  marg = list(
    l = 80,
    r = 20,
    t = 100,
    b = 60
  )
  
  mark = list(
    color = "lightblue",
    size = 10,
    line = list(
      color = "grey",
      width = 2
    )
  )
  
  fig = plot_ly(calib, x = ~ln_DTp, y = ~ln_CCSp, 
                name = 'data points', 
                type = 'scatter', 
                mode = "markers",
                marker = mark,
                text = ~paste(' m/z: ', round(calib$mz,4), '\n CCS (lit.):', calib$CCSN2,'Ų \n CCS (mes.):', calib$CCS_mes, 'Ų \n Difference:', calib$pct_diff_th, '%')) %>%
    add_trace(y = ~ln_CCS_th,
              name = 'calibration',
              mode = "lines+markers",
              marker = list(color = rgb(0.5,0.5,0.5,1), opacity=0),
              text = F,
              color = I("black")) %>%
    layout(title = titre, 
           xaxis = axis_x, 
           yaxis = axis_y, 
           margin = marg, 
           plot_bgcolor = rgb(225, 225, 225, maxColorValue = 255),
           showlegend = FALSE)
  
  
  # Apply calibration
  
  for(i in 1:length((data_as_list[[1]]))){
    
    # calculating the CCS
    
    drift_as_list[[2]][[i]][,2] = drift_as_list[[2]][[i]][,2]-0.001*delay*drift_as_list[[2]][[i]][,1]^0.5
    drift_as_list[[2]][[i]][,2] = drift_as_list[[2]][[i]][,2]^slope*drift_as_list[[2]][[i]][,4]*(1/(drift_as_list[[2]][[i]][,1]-(drift_as_list[[2]][[i]][,4]*cat_mass))+1/gas_mass)^0.5
    drift_as_list[[2]][[i]][,2] = (drift_as_list[[2]][[i]][,2]*coef_exp)/100
    
  }
  group = list(drift_as_list, fig)
  return(group)
}

#######################################################################################
write_mgf_IM = function(CCS_as_list, metadata, file_name, path_out){
  
  # Set the file name before writing
  output_file=paste(path_out, file_name, "_IM.mgf", sep="")
  
  # increment the pseudo mgf
  for(i in 1:length(CCS_as_list[[1]])){
    
    # browsing metadata
    ID = which(metadata$file_name==CCS_as_list[[1]][[i]])
    p_mass = metadata$mass[ID]
    z = metadata$charge[ID]
    
    # Writing the *.mgf
    text = "BEGIN IONS"
    text = append(text,paste("PEPMASS", p_mass, sep="="))
    text = append(text,paste("CHARGE",z,sep="="))
    text = append(text,"MSLEVEL=2")
    text = append(text,paste("SAMPLE",CCS_as_list[[1]][[i]],sep="="))
    text = append(text,paste("Title: Scan#: ", i, ", RT: 0.01 min", sep=""))
    for(j in 1:nrow(CCS_as_list[[2]][[i]])){
      line = paste((as.numeric(CCS_as_list[[2]][[i]][j,2])), 
                   CCS_as_list[[2]][[i]][j,3], sep=" ")
      text = append(text, line)
    }
    text = append(text, "END IONS")
    text = append(text, ' ')
    
    # Writing the file itself
    df_text=data.frame(text)
    write_tsv(df_text, output_file,  append = TRUE, col_names = FALSE)
  }
}

#######################################################################################
write_mgf_MZ = function(data_as_list, metadata, file_name, path_out){
  
  # Set the file name before writing
  output_file=paste(path_out, file_name, "_MS.mgf", sep="")
  
  # increment the pseudo mgf
  for(i in 1:length(data_as_list[[1]])){
    
    # browsing metadata
    ID = which(metadata$file_name==data_as_list[[1]][[i]])
    p_mass = metadata$mass[ID]
    z = metadata$charge[ID]
    
    # recreating the mz mass spectrum (bin fragments, sum intensities of same mass at 0.1 Da, mean of mz values)
    test=data_as_list[[2]][[i]]
    test$bin=round(test[,1],1)
    test2=unique(test$bin)
    test2$into=sum(test[which(test$bin==test2$`row m/z`),3])
    
    for(k in 1:nrow(test2)){
      a = which(test$bin==test2$`row m/z`[k])
      test2$into[k]=sum(test[a,3])
      test2$mz[k]=sum(test[a,1])/length(a)
    }    
    
    # Writing the *.mgf
    text = "BEGIN IONS"
    text = append(text,paste("PEPMASS", p_mass, sep="="))
    text = append(text,paste("CHARGE",z,sep="="))
    text = append(text,"MSLEVEL=2")
    text = append(text,paste("SAMPLE",data_as_list[[1]][[i]],sep="="))
    text = append(text,paste("Title: Scan#: ", i, ", RT: 0.01 min", sep=""))
    for(j in 1:nrow(test2)){
      line = paste(round(test2$mz[j],4), test2$into[j], sep=" ")
      text = append(text, line)
    }
    text = append(text, "END IONS")
    text = append(text, ' ')
    
    # Writing the file itself
    df_text=data.frame(text)
    write_tsv(df_text, output_file,  append = TRUE, col_names = FALSE)
    
    rm(test)
    rm(test2)  
  }
}

#######################################################################################
write_csv_clusterindex = function(data_as_list, file_name, path_out){
  
  # Set the file name before writing
  output_file=paste(path_out, file_name, "_CI.csv", sep="")
  
  # create data frame containing CI and sample names
  df_CI = matrix(data_as_list[[1]]) %>% unlist() %>% data.frame()
  colnames(df_CI) = c("file_name")
  df_CI$CI = 1:nrow(df_CI)
  
  write_csv(df_CI, output_file) 
}
