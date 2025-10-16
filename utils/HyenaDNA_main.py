import torch
import os,sys

import transformers
from transformers import PreTrainedModel, AutoModelForCausalLM, PretrainedConfig

sys.path.append('./HyenaDNA')
from huggingface import *
from standalone_hyenadna import *

max_length=16386
device='cpu'

# load pretrained model
model = HyenaDNAPreTrainedModel.from_pretrained(
    path="./HyenaDNA/",
    model_name='16k_enhancers_cohn',
    use_head=True,n_classes=2
)

# This two layer's parameter can not be load correctly automatically
checkpoint=torch.load("./HyenaDNA/16k_enhancers_cohn/weights.ckpt", map_location="cpu")
model.state_dict()['head.output_transform.weight'].copy_(checkpoint['state_dict']['decoder.0.output_transform.weight'])
model.state_dict()['head.output_transform.bias'].copy_(checkpoint['state_dict']['decoder.0.output_transform.bias'])

# create tokenizer, no training involved :)
tokenizer = CharacterTokenizer(
    characters=['A', 'C', 'G', 'T', 'N'],  # add DNA characters
    return_tensor='pt',
    model_max_length=max_length,
    padding='max_length',
    truncation=True
)

def hyena_inference(sequence,batch_size=1,device='cpu'):
    
    model.to(device)
    model.eval()
    
    with torch.inference_mode():
        if isinstance(sequence,str):
            tok_seq = tokenizer(sequence)["input_ids"]
            tok_seq = torch.LongTensor(tok_seq).unsqueeze(0).to(device)
            prob = model(tok_seq)
        elif isinstance(sequence,list):
            prob = []
            tok_seq = torch.concatenate([ torch.LongTensor(tokenizer(s)["input_ids"]).unsqueeze(0) for s in sequence ],dim=0).to(device)
            for i in range(0,tok_seq.shape[0],batch_size):
                if i+batch_size>=tok_seq.shape[0]:
                    prob.append(model(tok_seq[i:,:]))
                else:
                    prob.append(model(tok_seq[i:i+batch_size,:]))
            prob = torch.vstack(prob)
    
    return (prob[:,0]-prob[:,1]+1).cpu().numpy()/2*10