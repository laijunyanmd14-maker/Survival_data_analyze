#資料讀取
data = read.csv("melanoma.csv",header = TRUE)

###資料前處理
#設定censoring indicator
for(i in 1:nrow(data)){
  if(data[i,"Indicator"]==1){
    data[i,"delta"] = 1
  }else {
    data[i,"delta"] = 0
  }
}

#設定兩個dummy variables表示腫瘤厚度範圍
data$skin.1 = 0
data$skin.2 = 0

for(i in 1:nrow(data)){
  if(data[i,"Tumor.Thickness"]<=5 & data[i,"Tumor.Thickness"]>2){
    data[i,"skin.1"] = 1
  }else if(data[i,"Tumor.Thickness"]>5){
    data[i,"skin.2"] = 1
  }
}

#根據性別將資料分開
male = data[data$Sex==1,]#男性
female = data[data$Sex==0,]#女性

###敘述統計
#樣本數
table(data$Sex)

#定義比較函數(使用t檢定)
compare_by_sex_with_t_test = function(data,name){
  
  male = data[data$Sex==1,]
  female = data[data$Sex==0,]
  
  male_mean = mean(male[,name])
  male_sd = sd(male[,name])
  female_mean = mean(female[,name])
  female_sd = sd(female[,name])
  t_result = t.test(male[,name],female[,name],equal = FALSE)
  
  result = list(
    varible = name,
    male_mean = male_mean,
    male_sd = male_sd,
    female_mean = female_mean,
    female_sd = female_sd,
    t_result = t_result
    
  )
  
  return(result)
  
}

#追蹤時間分布與檢定
compare_by_sex_with_t_test(data,"tilde.T")

#年齡分布與檢定
compare_by_sex_with_t_test(data,"Age")

#比較腫瘤厚度比率
male_2 = nrow(male[male$Tumor.Thickness<=2,])
female_2 = nrow(female[female$Tumor.Thickness<=2,])
prop.test(x = c(male_2,female_2),n = c(nrow(male),nrow(female)))#男女腫瘤小於2之比率

male_3 = nrow(male[male$Tumor.Thickness>2 & male$Tumor.Thickness<=5,])
female_3 = nrow(female[female$Tumor.Thickness>2 & female$Tumor.Thickness<=5,])
prop.test(x = c(male_3,female_3),n = c(nrow(male),nrow(female)))#男女腫瘤介於2到5之比率

male_5 = nrow(male[male$Tumor.Thickness>5,])
female_5 = nrow(female[female$Tumor.Thickness>5,])
prop.test(x = c(male_5,female_5),n = c(nrow(male),nrow(female)))#男女腫瘤大於5之比率

#比較潰瘍比率
male_ulceration = nrow(male[male$Ulceration==1,])
female_ulceration = nrow(female[female$Ulceration==1,])
prop.test(x = c(male_ulceration,female_ulceration),n = c(nrow(male),nrow(female)))

#死亡原因比率
prop.table(table(male$Indicator))#男性
prop.table(table(female$Indicator))#女性

##畫KM並分析結果
# 載入套件
library(survival)
library(survminer)

# 建立生存物件
surv_obj <- Surv(time = data$tilde.T, event = data$delta)

# Kaplan–Meier，依性別分層
fit <- survfit(surv_obj ~ Sex, data = data)

# 繪製 KM 曲線
ggsurvplot(
  fit,
  data = data,
  pval = TRUE,            # 顯示 log-rank 檢定 p 值
  conf.int = TRUE,        # 顯示信賴區間
  risk.table = TRUE,      # 顯示風險人數表
  legend.labs = c("Female", "Male"),
  legend.title = "Sex",
  xlab = "Time",
  ylab = "Survival probability"
)

##檢驗生存曲線的差異
# log-rank test
survdiff(Surv(tilde.T, delta) ~ Sex, data = data)

## Kaplan–Meier, 分組: Sex + Ulceration
fit2 <- survfit(surv_obj ~ Sex + Ulceration, data = data)

# 繪製 KM 曲線
ggsurvplot(
  fit2,
  data = data,
  pval = TRUE,           # 顯示 log-rank p-value
  conf.int = TRUE,       # 顯示信賴區間
  risk.table = TRUE,     # 顯示風險人數表
  legend.title = "Sex + Ulceration",
  legend.labs = c("Female, No Ulcer", "Female, Ulcer",
                  "Male, No Ulcer", "Male, Ulcer"),
  xlab = "Time",
  ylab = "Survival probability",
  lwd = 1
)

##檢定性別與潰瘍與否對於生存曲線的影響
# 將性別與潰瘍與否合併成一個變數
data$group <- interaction(data$Sex, data$Ulceration)

# 檢定四組生存差異
survdiff(Surv(tilde.T, delta) ~ group, data = data)

##用Weibull加速失效時間模型分析黑色素瘤的死亡時間
# 建立 Weibull AFT 模型
aft_model <- survreg(
  surv_obj ~ Sex + Age + Ulceration + skin.1 + skin.2 + Sex:Ulceration,
  data = data,
  dist = "weibull"
)

summary(aft_model)

#腫瘤厚度的效果
exp(coef(aft_model)["skin.1"])
exp(coef(aft_model)["skin.2"])

#潰瘍的效果
exp(coef(aft_model)["Ulceration"])

#性別與潰瘍的交互作用效果
#建立無交互作用的模型
aft_model_no_int <- survreg(
  surv_obj ~ Sex + Age + Ulceration + skin.1 + skin.2,
  data = data,
  dist = "weibull"
)

#用 Likelihood Ratio Test 比較有無交互作用的模型
anova(aft_model_no_int, aft_model, test="Chisq")

##不考慮性別與潰瘍的交互作用而建立的模型
summary(aft_model_no_int)

#腫瘤厚度的效果
exp(coef(aft_model_no_int)["skin.1"])
exp(coef(aft_model_no_int)["skin.2"])

##計算男性、60歲、有潰瘍、腫瘤厚度3.5毫米的5年存活率和累積風險
# 取 coef 與 scale
coef <- coef(aft_model_no_int)
scale <- aft_model_no_int$scale  # σ

# X 值
Sex = 1 #男性
Age = 60 #60歲
Ulceration = 1 #有潰瘍
skin.1 = 1 #腫瘤厚度介於2-5mm
skin.2 = 0

# 線性預測值
eta <- coef["(Intercept)"] +
  coef["Sex"]*Sex +
  coef["Age"]*Age +
  coef["Ulceration"]*Ulceration +
  coef["skin.1"]*skin.1 +
  coef["skin.2"]*skin.2

# λ = exp(η)
lambda <- exp(eta)

# k = 1/σ
k <- 1/scale

time <- 5*365  # 5 年換成天數
lambda <- exp(eta)

S_5yr <- exp(-(time/lambda)^k)  # 5年存活率
H_5yr <- -log(S_5yr)         # 累積風險

eta <- as.numeric(eta)
S_5yr <- as.numeric(S_5yr)
H_5yr <- as.numeric(H_5yr)

S_5yr
H_5yr