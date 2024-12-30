import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score

fig, axs = plt.subplots(1, 2, figsize=(14, 7))

colors = ['#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#999999']

labels = [
    'proposed model', 'RPI-GGCN', 'LPICGAE', 'RPITER', 
    'GATLGEMF', 'LncPro', 'RPI-SE', 'LPI-Pred', 'RPI-CapsuleGAN'
]
# Here first the other methods are retrained and got their results. Then we named the results of the methods in .csv form as file 1-file 9 in order to make it suitable to read them in for loop.
for i in range(9):

    df = pd.read_csv(f'file{i + 1}.csv')

    true_labels = df['label']
    scores = df['pred']
    
    fpr, tpr, _ = roc_curve(true_labels, scores)
    roc_auc = auc(fpr, tpr)  
    axs[0].plot(fpr, tpr, color=colors[i], lw=1, label=f'{labels[i]} (AUC = {roc_auc:.2f})')

    precision, recall, _ = precision_recall_curve(true_labels, scores)
    auprc = average_precision_score(true_labels, scores) 
    axs[1].plot(recall, precision, color=colors[i], lw=1, label=f'{labels[i]} (AUPRC = {auprc:.2f})')

#axs[0].plot([0, 1], [0, 1], color='navy', lw=1, linestyle='--')
axs[0].set_xlim([0.0, 1.0])
axs[0].set_ylim([0.0, 1.05])
axs[0].set_xlabel('False Positive Rate')
axs[0].set_ylabel('True Positive Rate')
axs[0].set_title('Receiver Operating Characteristic (ROC) Curve')
axs[0].legend(loc='lower right')

axs[1].set_xlim([0.0, 1.0])
axs[1].set_ylim([0.0, 1.05])
axs[1].set_xlabel('Recall')
axs[1].set_ylabel('Precision')
axs[1].set_title('Precision-Recall Curve')
axs[1].legend(loc='lower left')

plt.tight_layout()
plt.show()
