from torch import nn

class SimpleMLP(nn.Module):
    def __init__(self, input_size=2560, output_size=1):
        super(SimpleMLP, self).__init__()

        self.project = nn.Sequential(
                nn.Linear(input_size, input_size),
                nn.ReLU(),
                nn.Dropout(0.2),
                nn.Linear(input_size, output_size)
            )
        
    def forward(self, x):
        return self.project(x)