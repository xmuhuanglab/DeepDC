import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.utils.data import TensorDataset, DataLoader
from tqdm import tqdm
import numpy as np

base_dict = {'A':[1,0,0,0],'C':[0,1,0,0],'G':[0,0,1,0],'T':[0,0,0,1],'N':[0,0,0,0]}

def one_hot(x):return np.array([base_dict[b] for b in x])

def mish(x):return x * torch.tanh(F.softplus(x))

class DistillatedCRISPRedict(nn.Module):
    def __init__(self):
        super().__init__()
        self.fc0 = nn.Linear(120, 1024)
        self.norm0 = nn.LayerNorm(1024)
        self.fc1 = nn.Linear(1024, 1024)
        self.norm1 = nn.LayerNorm(1024)
        self.fc2 = nn.Linear(1024, 512)
        self.norm2 = nn.LayerNorm(512)
        self.fc3 = nn.Linear(512, 256)
        self.norm3 = nn.LayerNorm(256)
        self.output = nn.Linear(256, 1)
        
    def forward(self, x):
        x = x.view(x.size(0),-1)
        
        x = self.fc0(x)
        x = self.norm0(x)
        x = mish(x)
        
        x = self.fc1(x)
        x = self.norm1(x)
        x = mish(x)
        
        x = self.fc2(x)
        x = self.norm2(x)
        x = mish(x)
        
        x = self.fc3(x)
        x = self.norm3(x)
        x = mish(x)        
        
        return self.output(x).flatten()
    
    def fit(self, x, y, epochs=20, batch_size=256, lr=1e-5, device='cpu'):
        if not isinstance(x, torch.Tensor):
            x = torch.FloatTensor(x).to(device)
        if not isinstance(y, torch.Tensor):
            y = torch.FloatTensor(y).to(device)

        dataset = TensorDataset(x, y)
        dataloader = DataLoader(dataset, batch_size=batch_size, shuffle=True)

        optimizer = torch.optim.Adam(self.parameters(), lr=lr)
        criterion = nn.MSELoss() 

        self.train()
        self.to(device)
        for epoch in range(epochs):
            epoch_loss = 0.0
            for batch_x, batch_y in tqdm(dataloader):
                optimizer.zero_grad()

                outputs = self(batch_x)
                loss = criterion(outputs, batch_y)

                loss.backward()
                optimizer.step()
                
                epoch_loss += loss.item()

            print(f'Epoch [{epoch+1}/{epochs}], Loss: {epoch_loss/len(dataloader):.4f}')
            
    def predict(self, x, batch_size=256):
        if not isinstance(x, torch.Tensor):
            x = torch.FloatTensor(x)
        dataset = TensorDataset(x)
        dataloader = DataLoader(dataset, batch_size=batch_size)
        all_preds = []
        self.eval()
        with torch.no_grad():
            for batch_x in dataloader:
                preds = self(batch_x[0])
                all_preds.append(preds)
        return torch.cat(all_preds).numpy()
    
def predict(sequences,model_path):
    model = DistillatedCRISPRedict()
    model.load_state_dict(torch.load(model_path,map_location=torch.device('cpu')))
    x = np.stack( [ one_hot(seq) for seq in sequences] )
    y_pred = model.predict(x)
    return y_pred
