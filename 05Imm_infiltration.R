# This code loads, reshape and filter immune infiltration preprocesses data
# and creates a box plot to visualize the cell composition in the tumor microenvironment.
#
# 1. Data Loading and Preparation:
#
# Read the "estimation_matrix.csv" file into the 'imm' data frame, transpose it, and clean up.
# Read the "mut_df.csv" file into the 'mut_df' data frame and remove the first column.
# 
# 2. Data Manipulation and Merging:
#   
# Extract the 'Amino_Acid_Change' column from 'mut_df' into 'group_list'.
# PJoin between 'imm' and 'mut_df' based on 'ID' and 'X' columns, creating 'TME_data'.
# 
# 3. Data Filtering:
#   
# Filter out unwanted rows with specific 'Celltype' values from the reshaped 'TME_data'.
# 
# 4. Plot Customization and Creation:
#   
# Load the 'ggpubr' library for creating publication-ready plots.
# Define a custom theme for the upcoming plot.
# 
# 5. Box Plot Generation:
#   
# Create a box plot using 'ggplot' with data from the reshaped 'TME_data'.
# Label the plot axes and title.
# Customize the appearance of the box plot, including grouping by 'Group' and coloring.
# Apply the custom theme and add significance labels using 'stat_compare_means'.
# 
# 

# Set working directory to the specified path
# setwd("./data/")

# Read and Transposethe file into the 'imm' data frame
imm <- read.csv("./estimation_matrix.csv",header = T)
imm <- as.data.frame(t(imm))
colnames(imm) <- imm[1,]# Use the first row of 'imm' as column names
imm <- imm[-1,]
imm$ID <- rownames(imm)# Create a new column 'ID' in 'imm' containing row names

# Read the file "mut_df.csv" with column data from 'data_processing.R' 
mut_df <- read.csv("./data/02_exps/TCGA/mut_df.csv",header = T,sep = ',')
mut_df <- mut_df[,-1]
group_list <- mut_df$Amino_Acid_Change # Extract the 'Amino_Acid_Change' column into 'group_list'
db <- left_join(imm,mut_df,by=c('ID'='X')) # Left join 'imm' and 'mut_df' data frames using the 'ID' and 'X' columns
TME_data <- db


library(dplyr)
library(reshape2)
# Rename columns in 'TME_data' data frame
names(TME_data)[names(TME_data) == "Amino_Acid_Change"] <- "group"
names(TME_data)[names(TME_data) == "ID"] <- "sample"

# Reshape 'TME_data' from wide to long format
TME_New = melt(TME_data,measure.vars = colnames(TME_data),id.vars = c("group","sample"))
# Rename the columns of 'TME_New'
colnames(TME_New)=c("Group","Sample","Celltype","Composition")
head(TME_New)

# plot order by median
plot_order = TME_New[TME_New$Group=="Del19",] %>% 
  group_by(Celltype) %>% 
  summarise(m = median(Composition)) %>% 
  arrange(desc(m)) %>% 
  pull(Celltype)

# Filter out unwanted rows with specific 'Celltype' values

TME_New$Celltype = factor(TME_New$Celltype,levels = plot_order)
TME_New$Composition <- as.numeric(TME_New$Composition)
print(class(TME_New$Composition))
TME_New$Composition <- round(TME_New$Composition,4)
TME_New <- TME_New[TME_New$Celltype != c("sample"), ]
TME_New <- TME_New[TME_New$Celltype != c("Sample_ID"), ]
TME_New <- TME_New[TME_New$Celltype != c("group"), ]


# install.packages("ggpubr")
library(ggpubr)

# Check if the condition 'T' is true (which it typically is)
if(T){
  # Define custom theme for the plot
  mytheme <- theme(plot.title = element_text(size = 12,color="black",hjust = 0.5),
                   axis.title = element_text(size = 12,color ="black"), 
                   axis.text = element_text(size= 12,color = "black"),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.minor.x = element_blank(),
                   axis.text.x = element_text(angle = 45, hjust = 1 ),
                   panel.grid=element_blank(),
                   legend.position = "top",
                   legend.text = element_text(size= 12),
                   legend.title= element_text(size= 12)
  ) 
}

# Create a box plot using 'ggplot' with data from 'TME_New'
box_TME <- ggplot(TME_New, aes(x = Celltype, y = Composition))+ 
  labs(y="Cell composition",x= NULL,title = "TME Cell composition")+  
  geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
  scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
  theme_classic() + mytheme + 
  stat_compare_means(aes(group =  Group),
                     label = "p.signif",
                     method = "wilcox.test",
                     hide.ns = T)

# Display the box plot
box_TME

# Save the box plot as a PDF file
ggsave("./TCGA_HNSCC_TME.pdf",box_TME,height=15,width=35,unit="cm")











