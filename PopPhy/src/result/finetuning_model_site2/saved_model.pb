��
��
^
AssignVariableOp
resource
value"dtype"
dtypetype"
validate_shapebool( �
�
BiasAdd

value"T	
bias"T
output"T""
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
8
Const
output"dtype"
valuetensor"
dtypetype
�
Conv2D

input"T
filter"T
output"T"
Ttype:	
2"
strides	list(int)"
use_cudnn_on_gpubool(",
paddingstring:
SAMEVALIDEXPLICIT""
explicit_paddings	list(int)
 "-
data_formatstringNHWC:
NHWCNCHW" 
	dilations	list(int)

$
DisableCopyOnRead
resource�
.
Identity

input"T
output"T"	
Ttype
u
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:
2	
�
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool("
allow_missing_filesbool( �

NoOp
M
Pack
values"T*N
output"T"
Nint(0"	
Ttype"
axisint 
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetype�
E
Relu
features"T
activations"T"
Ttype:
2	
[
Reshape
tensor"T
shape"Tshape
output"T"	
Ttype"
Tshapetype0:
2	
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
?
Select
	condition

t"T
e"T
output"T"	
Ttype
H
ShardedFilename
basename	
shard

num_shards
filename
�
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring ��
@
StaticRegexFullMatch	
input

output
"
patternstring
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 
�
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 �"serve*2.12.02v2.12.0-rc1-12-g0db597d0d758��
^
countVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namecount
W
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
dtype0
^
totalVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nametotal
W
total/Read/ReadVariableOpReadVariableOptotal*
_output_shapes
: *
dtype0
b
count_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_1
[
count_1/Read/ReadVariableOpReadVariableOpcount_1*
_output_shapes
: *
dtype0
b
total_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_1
[
total_1/Read/ReadVariableOpReadVariableOptotal_1*
_output_shapes
: *
dtype0
~
Adam/v/dense_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*$
shared_nameAdam/v/dense_3/bias
w
'Adam/v/dense_3/bias/Read/ReadVariableOpReadVariableOpAdam/v/dense_3/bias*
_output_shapes
:
*
dtype0
~
Adam/m/dense_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*$
shared_nameAdam/m/dense_3/bias
w
'Adam/m/dense_3/bias/Read/ReadVariableOpReadVariableOpAdam/m/dense_3/bias*
_output_shapes
:
*
dtype0
�
Adam/v/dense_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
: 
*&
shared_nameAdam/v/dense_3/kernel

)Adam/v/dense_3/kernel/Read/ReadVariableOpReadVariableOpAdam/v/dense_3/kernel*
_output_shapes

: 
*
dtype0
�
Adam/m/dense_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
: 
*&
shared_nameAdam/m/dense_3/kernel

)Adam/m/dense_3/kernel/Read/ReadVariableOpReadVariableOpAdam/m/dense_3/kernel*
_output_shapes

: 
*
dtype0
x
Adam/v/fc_0/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *!
shared_nameAdam/v/fc_0/bias
q
$Adam/v/fc_0/bias/Read/ReadVariableOpReadVariableOpAdam/v/fc_0/bias*
_output_shapes
: *
dtype0
x
Adam/m/fc_0/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *!
shared_nameAdam/m/fc_0/bias
q
$Adam/m/fc_0/bias/Read/ReadVariableOpReadVariableOpAdam/m/fc_0/bias*
_output_shapes
: *
dtype0
�
Adam/v/fc_0/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�. *#
shared_nameAdam/v/fc_0/kernel
z
&Adam/v/fc_0/kernel/Read/ReadVariableOpReadVariableOpAdam/v/fc_0/kernel*
_output_shapes
:	�. *
dtype0
�
Adam/m/fc_0/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�. *#
shared_nameAdam/m/fc_0/kernel
z
&Adam/m/fc_0/kernel/Read/ReadVariableOpReadVariableOpAdam/m/fc_0/kernel*
_output_shapes
:	�. *
dtype0
|
Adam/v/conv_1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/v/conv_1/bias
u
&Adam/v/conv_1/bias/Read/ReadVariableOpReadVariableOpAdam/v/conv_1/bias*
_output_shapes
: *
dtype0
|
Adam/m/conv_1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/m/conv_1/bias
u
&Adam/m/conv_1/bias/Read/ReadVariableOpReadVariableOpAdam/m/conv_1/bias*
_output_shapes
: *
dtype0
�
Adam/v/conv_1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:  *%
shared_nameAdam/v/conv_1/kernel
�
(Adam/v/conv_1/kernel/Read/ReadVariableOpReadVariableOpAdam/v/conv_1/kernel*&
_output_shapes
:  *
dtype0
�
Adam/m/conv_1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:  *%
shared_nameAdam/m/conv_1/kernel
�
(Adam/m/conv_1/kernel/Read/ReadVariableOpReadVariableOpAdam/m/conv_1/kernel*&
_output_shapes
:  *
dtype0
|
Adam/v/conv_0/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/v/conv_0/bias
u
&Adam/v/conv_0/bias/Read/ReadVariableOpReadVariableOpAdam/v/conv_0/bias*
_output_shapes
: *
dtype0
|
Adam/m/conv_0/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *#
shared_nameAdam/m/conv_0/bias
u
&Adam/m/conv_0/bias/Read/ReadVariableOpReadVariableOpAdam/m/conv_0/bias*
_output_shapes
: *
dtype0
�
Adam/v/conv_0/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape: *%
shared_nameAdam/v/conv_0/kernel
�
(Adam/v/conv_0/kernel/Read/ReadVariableOpReadVariableOpAdam/v/conv_0/kernel*&
_output_shapes
: *
dtype0
�
Adam/m/conv_0/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape: *%
shared_nameAdam/m/conv_0/kernel
�
(Adam/m/conv_0/kernel/Read/ReadVariableOpReadVariableOpAdam/m/conv_0/kernel*&
_output_shapes
: *
dtype0
n
learning_rateVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namelearning_rate
g
!learning_rate/Read/ReadVariableOpReadVariableOplearning_rate*
_output_shapes
: *
dtype0
f
	iterationVarHandleOp*
_output_shapes
: *
dtype0	*
shape: *
shared_name	iteration
_
iteration/Read/ReadVariableOpReadVariableOp	iteration*
_output_shapes
: *
dtype0	
j
	fc_0/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	fc_0/bias
c
fc_0/bias/Read/ReadVariableOpReadVariableOp	fc_0/bias*
_output_shapes
: *
dtype0
s
fc_0/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	�. *
shared_namefc_0/kernel
l
fc_0/kernel/Read/ReadVariableOpReadVariableOpfc_0/kernel*
_output_shapes
:	�. *
dtype0
n
conv_1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameconv_1/bias
g
conv_1/bias/Read/ReadVariableOpReadVariableOpconv_1/bias*
_output_shapes
: *
dtype0
~
conv_1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:  *
shared_nameconv_1/kernel
w
!conv_1/kernel/Read/ReadVariableOpReadVariableOpconv_1/kernel*&
_output_shapes
:  *
dtype0
n
conv_0/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameconv_0/bias
g
conv_0/bias/Read/ReadVariableOpReadVariableOpconv_0/bias*
_output_shapes
: *
dtype0
~
conv_0/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameconv_0/kernel
w
!conv_0/kernel/Read/ReadVariableOpReadVariableOpconv_0/kernel*&
_output_shapes
: *
dtype0
p
dense_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:
*
shared_namedense_3/bias
i
 dense_3/bias/Read/ReadVariableOpReadVariableOpdense_3/bias*
_output_shapes
:
*
dtype0
x
dense_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
: 
*
shared_namedense_3/kernel
q
"dense_3/kernel/Read/ReadVariableOpReadVariableOpdense_3/kernel*
_output_shapes

: 
*
dtype0
�
serving_default_input_5Placeholder*/
_output_shapes
:���������
'*
dtype0*$
shape:���������
'
�
StatefulPartitionedCallStatefulPartitionedCallserving_default_input_5conv_0/kernelconv_0/biasconv_1/kernelconv_1/biasfc_0/kernel	fc_0/biasdense_3/kerneldense_3/bias*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *,
f'R%
#__inference_signature_wrapper_11169

NoOpNoOp
�N
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*�N
value�NB�N B�N
�
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*	&call_and_return_all_conditional_losses

_default_save_signature
	optimizer

signatures*
�
layer-0
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses* 
�
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer-3
layer_with_weights-2
layer-4
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses*
�
	variables
 trainable_variables
!regularization_losses
"	keras_api
#__call__
*$&call_and_return_all_conditional_losses

%kernel
&bias*
<
'0
(1
)2
*3
+4
,5
%6
&7*
<
'0
(1
)2
*3
+4
,5
%6
&7*
* 
�
-non_trainable_variables

.layers
/metrics
0layer_regularization_losses
1layer_metrics
	variables
trainable_variables
regularization_losses
__call__

_default_save_signature
*	&call_and_return_all_conditional_losses
&	"call_and_return_conditional_losses*
6
2trace_0
3trace_1
4trace_2
5trace_3* 
6
6trace_0
7trace_1
8trace_2
9trace_3* 
* 
�
:
_variables
;_iterations
<_learning_rate
=_index_dict
>
_momentums
?_velocities
@_update_step_xla*

Aserving_default* 
�
B	variables
Ctrainable_variables
Dregularization_losses
E	keras_api
F__call__
*G&call_and_return_all_conditional_losses* 
* 
* 
* 
�
Hnon_trainable_variables

Ilayers
Jmetrics
Klayer_regularization_losses
Llayer_metrics
	variables
trainable_variables
regularization_losses
__call__
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses* 
6
Mtrace_0
Ntrace_1
Otrace_2
Ptrace_3* 
6
Qtrace_0
Rtrace_1
Strace_2
Ttrace_3* 
�
U	variables
Vtrainable_variables
Wregularization_losses
X	keras_api
Y__call__
*Z&call_and_return_all_conditional_losses
[_random_generator* 
�
\	variables
]trainable_variables
^regularization_losses
_	keras_api
`__call__
*a&call_and_return_all_conditional_losses

'kernel
(bias
 b_jit_compiled_convolution_op*
�
c	variables
dtrainable_variables
eregularization_losses
f	keras_api
g__call__
*h&call_and_return_all_conditional_losses

)kernel
*bias
 i_jit_compiled_convolution_op*
�
j	variables
ktrainable_variables
lregularization_losses
m	keras_api
n__call__
*o&call_and_return_all_conditional_losses* 
�
p	variables
qtrainable_variables
rregularization_losses
s	keras_api
t__call__
*u&call_and_return_all_conditional_losses

+kernel
,bias*
.
'0
(1
)2
*3
+4
,5*
.
'0
(1
)2
*3
+4
,5*
,
v0
w1
x2
y3
z4
{5* 
�
|non_trainable_variables

}layers
~metrics
layer_regularization_losses
�layer_metrics
	variables
trainable_variables
regularization_losses
__call__
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses*
:
�trace_0
�trace_1
�trace_2
�trace_3* 
:
�trace_0
�trace_1
�trace_2
�trace_3* 

%0
&1*

%0
&1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
	variables
 trainable_variables
!regularization_losses
#__call__
*$&call_and_return_all_conditional_losses
&$"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
^X
VARIABLE_VALUEdense_3/kernel6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUEdense_3/bias4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE*
MG
VARIABLE_VALUEconv_0/kernel&variables/0/.ATTRIBUTES/VARIABLE_VALUE*
KE
VARIABLE_VALUEconv_0/bias&variables/1/.ATTRIBUTES/VARIABLE_VALUE*
MG
VARIABLE_VALUEconv_1/kernel&variables/2/.ATTRIBUTES/VARIABLE_VALUE*
KE
VARIABLE_VALUEconv_1/bias&variables/3/.ATTRIBUTES/VARIABLE_VALUE*
KE
VARIABLE_VALUEfc_0/kernel&variables/4/.ATTRIBUTES/VARIABLE_VALUE*
IC
VARIABLE_VALUE	fc_0/bias&variables/5/.ATTRIBUTES/VARIABLE_VALUE*
* 

0
1
2*

�0
�1*
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
�
;0
�1
�2
�3
�4
�5
�6
�7
�8
�9
�10
�11
�12
�13
�14
�15
�16*
SM
VARIABLE_VALUE	iteration0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUElearning_rate3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
* 
D
�0
�1
�2
�3
�4
�5
�6
�7*
D
�0
�1
�2
�3
�4
�5
�6
�7*
* 
* 
* 
* 
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
B	variables
Ctrainable_variables
Dregularization_losses
F__call__
*G&call_and_return_all_conditional_losses
&G"call_and_return_conditional_losses* 

�trace_0
�trace_1* 

�trace_0
�trace_1* 
* 
	
0* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
U	variables
Vtrainable_variables
Wregularization_losses
Y__call__
*Z&call_and_return_all_conditional_losses
&Z"call_and_return_conditional_losses* 

�trace_0
�trace_1* 

�trace_0
�trace_1* 
* 

'0
(1*

'0
(1*

v0
w1* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
\	variables
]trainable_variables
^regularization_losses
`__call__
*a&call_and_return_all_conditional_losses
&a"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
* 

)0
*1*

)0
*1*

x0
y1* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
c	variables
dtrainable_variables
eregularization_losses
g__call__
*h&call_and_return_all_conditional_losses
&h"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
* 
* 
* 
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
j	variables
ktrainable_variables
lregularization_losses
n__call__
*o&call_and_return_all_conditional_losses
&o"call_and_return_conditional_losses* 

�trace_0* 

�trace_0* 

+0
,1*

+0
,1*

z0
{1* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
p	variables
qtrainable_variables
rregularization_losses
t__call__
*u&call_and_return_all_conditional_losses
&u"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 

�trace_0* 

�trace_0* 

�trace_0* 

�trace_0* 

�trace_0* 

�trace_0* 
* 
'
0
1
2
3
4*
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
<
�	variables
�	keras_api

�total

�count*
M
�	variables
�	keras_api

�total

�count
�
_fn_kwargs*
_Y
VARIABLE_VALUEAdam/m/conv_0/kernel1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/v/conv_0/kernel1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUE*
]W
VARIABLE_VALUEAdam/m/conv_0/bias1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUE*
]W
VARIABLE_VALUEAdam/v/conv_0/bias1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/m/conv_1/kernel1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/v/conv_1/kernel1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUE*
]W
VARIABLE_VALUEAdam/m/conv_1/bias1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUE*
]W
VARIABLE_VALUEAdam/v/conv_1/bias1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUE*
]W
VARIABLE_VALUEAdam/m/fc_0/kernel1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUE*
^X
VARIABLE_VALUEAdam/v/fc_0/kernel2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUE*
\V
VARIABLE_VALUEAdam/m/fc_0/bias2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUE*
\V
VARIABLE_VALUEAdam/v/fc_0/bias2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/m/dense_3/kernel2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUE*
a[
VARIABLE_VALUEAdam/v/dense_3/kernel2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/m/dense_3/bias2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUE*
_Y
VARIABLE_VALUEAdam/v/dense_3/bias2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUE*
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 

v0
w1* 
* 
* 
* 
* 
* 
* 

x0
y1* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 

z0
{1* 
* 
* 
* 
* 
* 
* 
* 
* 
* 

�0
�1*

�	variables*
UO
VARIABLE_VALUEtotal_14keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEcount_14keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE*

�0
�1*

�	variables*
SM
VARIABLE_VALUEtotal4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUE*
SM
VARIABLE_VALUEcount4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUE*
* 
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
�
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filenamedense_3/kerneldense_3/biasconv_0/kernelconv_0/biasconv_1/kernelconv_1/biasfc_0/kernel	fc_0/bias	iterationlearning_rateAdam/m/conv_0/kernelAdam/v/conv_0/kernelAdam/m/conv_0/biasAdam/v/conv_0/biasAdam/m/conv_1/kernelAdam/v/conv_1/kernelAdam/m/conv_1/biasAdam/v/conv_1/biasAdam/m/fc_0/kernelAdam/v/fc_0/kernelAdam/m/fc_0/biasAdam/v/fc_0/biasAdam/m/dense_3/kernelAdam/v/dense_3/kernelAdam/m/dense_3/biasAdam/v/dense_3/biastotal_1count_1totalcountConst*+
Tin$
"2 *
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *'
f"R 
__inference__traced_save_12060
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenamedense_3/kerneldense_3/biasconv_0/kernelconv_0/biasconv_1/kernelconv_1/biasfc_0/kernel	fc_0/bias	iterationlearning_rateAdam/m/conv_0/kernelAdam/v/conv_0/kernelAdam/m/conv_0/biasAdam/v/conv_0/biasAdam/m/conv_1/kernelAdam/v/conv_1/kernelAdam/m/conv_1/biasAdam/v/conv_1/biasAdam/m/fc_0/kernelAdam/v/fc_0/kernelAdam/m/fc_0/biasAdam/v/fc_0/biasAdam/m/dense_3/kernelAdam/v/dense_3/kernelAdam/m/dense_3/biasAdam/v/dense_3/biastotal_1count_1totalcount**
Tin#
!2*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� **
f%R#
!__inference__traced_restore_12160��
�O
�	
K__inference_finetuning_model_layer_call_and_return_conditional_losses_11391

inputsL
2sequential_2_conv_0_conv2d_readvariableop_resource: A
3sequential_2_conv_0_biasadd_readvariableop_resource: L
2sequential_2_conv_1_conv2d_readvariableop_resource:  A
3sequential_2_conv_1_biasadd_readvariableop_resource: C
0sequential_2_fc_0_matmul_readvariableop_resource:	�. ?
1sequential_2_fc_0_biasadd_readvariableop_resource: 8
&dense_3_matmul_readvariableop_resource: 
5
'dense_3_biasadd_readvariableop_resource:

identity��-conv_0/bias/Regularizer/L2Loss/ReadVariableOp�/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp�-conv_1/bias/Regularizer/L2Loss/ReadVariableOp�/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp�dense_3/BiasAdd/ReadVariableOp�dense_3/MatMul/ReadVariableOp�+fc_0/bias/Regularizer/L2Loss/ReadVariableOp�-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp�*sequential_2/conv_0/BiasAdd/ReadVariableOp�)sequential_2/conv_0/Conv2D/ReadVariableOp�*sequential_2/conv_1/BiasAdd/ReadVariableOp�)sequential_2/conv_1/Conv2D/ReadVariableOp�(sequential_2/fc_0/BiasAdd/ReadVariableOp�'sequential_2/fc_0/MatMul/ReadVariableOp�
)sequential_2/conv_0/Conv2D/ReadVariableOpReadVariableOp2sequential_2_conv_0_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype0�
sequential_2/conv_0/Conv2DConv2Dinputs1sequential_2/conv_0/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������# *
paddingVALID*
strides
�
*sequential_2/conv_0/BiasAdd/ReadVariableOpReadVariableOp3sequential_2_conv_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
sequential_2/conv_0/BiasAddBiasAdd#sequential_2/conv_0/Conv2D:output:02sequential_2/conv_0/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������# �
sequential_2/conv_0/ReluRelu$sequential_2/conv_0/BiasAdd:output:0*
T0*/
_output_shapes
:���������# �
)sequential_2/conv_1/Conv2D/ReadVariableOpReadVariableOp2sequential_2_conv_1_conv2d_readvariableop_resource*&
_output_shapes
:  *
dtype0�
sequential_2/conv_1/Conv2DConv2D&sequential_2/conv_0/Relu:activations:01sequential_2/conv_1/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:��������� *
paddingVALID*
strides
�
*sequential_2/conv_1/BiasAdd/ReadVariableOpReadVariableOp3sequential_2_conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
sequential_2/conv_1/BiasAddBiasAdd#sequential_2/conv_1/Conv2D:output:02sequential_2/conv_1/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:��������� �
sequential_2/conv_1/ReluRelu$sequential_2/conv_1/BiasAdd:output:0*
T0*/
_output_shapes
:��������� k
sequential_2/flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"����@  �
sequential_2/flatten/ReshapeReshape&sequential_2/conv_1/Relu:activations:0#sequential_2/flatten/Const:output:0*
T0*(
_output_shapes
:����������.�
'sequential_2/fc_0/MatMul/ReadVariableOpReadVariableOp0sequential_2_fc_0_matmul_readvariableop_resource*
_output_shapes
:	�. *
dtype0�
sequential_2/fc_0/MatMulMatMul%sequential_2/flatten/Reshape:output:0/sequential_2/fc_0/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� �
(sequential_2/fc_0/BiasAdd/ReadVariableOpReadVariableOp1sequential_2_fc_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
sequential_2/fc_0/BiasAddBiasAdd"sequential_2/fc_0/MatMul:product:00sequential_2/fc_0/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� t
sequential_2/fc_0/ReluRelu"sequential_2/fc_0/BiasAdd:output:0*
T0*'
_output_shapes
:��������� �
dense_3/MatMul/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource*
_output_shapes

: 
*
dtype0�
dense_3/MatMulMatMul$sequential_2/fc_0/Relu:activations:0%dense_3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
dense_3/BiasAdd/ReadVariableOpReadVariableOp'dense_3_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype0�
dense_3/BiasAddBiasAdddense_3/MatMul:product:0&dense_3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp2sequential_2_conv_0_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype0�
 conv_0/kernel/Regularizer/L2LossL2Loss7conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/kernel/Regularizer/mulMul(conv_0/kernel/Regularizer/mul/x:output:0)conv_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-conv_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp3sequential_2_conv_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_0/bias/Regularizer/L2LossL2Loss5conv_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/bias/Regularizer/mulMul&conv_0/bias/Regularizer/mul/x:output:0'conv_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp2sequential_2_conv_1_conv2d_readvariableop_resource*&
_output_shapes
:  *
dtype0�
 conv_1/kernel/Regularizer/L2LossL2Loss7conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/kernel/Regularizer/mulMul(conv_1/kernel/Regularizer/mul/x:output:0)conv_1/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-conv_1/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp3sequential_2_conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_1/bias/Regularizer/L2LossL2Loss5conv_1/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/bias/Regularizer/mulMul&conv_1/bias/Regularizer/mul/x:output:0'conv_1/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp0sequential_2_fc_0_matmul_readvariableop_resource*
_output_shapes
:	�. *
dtype0�
fc_0/kernel/Regularizer/L2LossL2Loss5fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
fc_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/kernel/Regularizer/mulMul&fc_0/kernel/Regularizer/mul/x:output:0'fc_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
+fc_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp1sequential_2_fc_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0|
fc_0/bias/Regularizer/L2LossL2Loss3fc_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: `
fc_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/bias/Regularizer/mulMul$fc_0/bias/Regularizer/mul/x:output:0%fc_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: g
IdentityIdentitydense_3/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������
�
NoOpNoOp.^conv_0/bias/Regularizer/L2Loss/ReadVariableOp0^conv_0/kernel/Regularizer/L2Loss/ReadVariableOp.^conv_1/bias/Regularizer/L2Loss/ReadVariableOp0^conv_1/kernel/Regularizer/L2Loss/ReadVariableOp^dense_3/BiasAdd/ReadVariableOp^dense_3/MatMul/ReadVariableOp,^fc_0/bias/Regularizer/L2Loss/ReadVariableOp.^fc_0/kernel/Regularizer/L2Loss/ReadVariableOp+^sequential_2/conv_0/BiasAdd/ReadVariableOp*^sequential_2/conv_0/Conv2D/ReadVariableOp+^sequential_2/conv_1/BiasAdd/ReadVariableOp*^sequential_2/conv_1/Conv2D/ReadVariableOp)^sequential_2/fc_0/BiasAdd/ReadVariableOp(^sequential_2/fc_0/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:���������
': : : : : : : : 2^
-conv_0/bias/Regularizer/L2Loss/ReadVariableOp-conv_0/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp2^
-conv_1/bias/Regularizer/L2Loss/ReadVariableOp-conv_1/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp2@
dense_3/BiasAdd/ReadVariableOpdense_3/BiasAdd/ReadVariableOp2>
dense_3/MatMul/ReadVariableOpdense_3/MatMul/ReadVariableOp2Z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOp+fc_0/bias/Regularizer/L2Loss/ReadVariableOp2^
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp2X
*sequential_2/conv_0/BiasAdd/ReadVariableOp*sequential_2/conv_0/BiasAdd/ReadVariableOp2V
)sequential_2/conv_0/Conv2D/ReadVariableOp)sequential_2/conv_0/Conv2D/ReadVariableOp2X
*sequential_2/conv_1/BiasAdd/ReadVariableOp*sequential_2/conv_1/BiasAdd/ReadVariableOp2V
)sequential_2/conv_1/Conv2D/ReadVariableOp)sequential_2/conv_1/Conv2D/ReadVariableOp2T
(sequential_2/fc_0/BiasAdd/ReadVariableOp(sequential_2/fc_0/BiasAdd/ReadVariableOp2R
'sequential_2/fc_0/MatMul/ReadVariableOp'sequential_2/fc_0/MatMul/ReadVariableOp:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�0
f
G__inference_sequential_3_layer_call_and_return_conditional_losses_11440

inputs
identity�_
random_color_affine_2/ShapeShapeinputs*
T0*
_output_shapes
::��s
)random_color_affine_2/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: u
+random_color_affine_2/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:u
+random_color_affine_2/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:�
#random_color_affine_2/strided_sliceStridedSlice$random_color_affine_2/Shape:output:02random_color_affine_2/strided_slice/stack:output:04random_color_affine_2/strided_slice/stack_1:output:04random_color_affine_2/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_maskn
,random_color_affine_2/random_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :n
,random_color_affine_2/random_uniform/shape/2Const*
_output_shapes
: *
dtype0*
value	B :n
,random_color_affine_2/random_uniform/shape/3Const*
_output_shapes
: *
dtype0*
value	B :�
*random_color_affine_2/random_uniform/shapePack,random_color_affine_2/strided_slice:output:05random_color_affine_2/random_uniform/shape/1:output:05random_color_affine_2/random_uniform/shape/2:output:05random_color_affine_2/random_uniform/shape/3:output:0*
N*
T0*
_output_shapes
:m
(random_color_affine_2/random_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *����m
(random_color_affine_2/random_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *���>�
2random_color_affine_2/random_uniform/RandomUniformRandomUniform3random_color_affine_2/random_uniform/shape:output:0*
T0*/
_output_shapes
:���������*
dtype0�
(random_color_affine_2/random_uniform/subSub1random_color_affine_2/random_uniform/max:output:01random_color_affine_2/random_uniform/min:output:0*
T0*
_output_shapes
: �
(random_color_affine_2/random_uniform/mulMul;random_color_affine_2/random_uniform/RandomUniform:output:0,random_color_affine_2/random_uniform/sub:z:0*
T0*/
_output_shapes
:����������
$random_color_affine_2/random_uniformAddV2,random_color_affine_2/random_uniform/mul:z:01random_color_affine_2/random_uniform/min:output:0*
T0*/
_output_shapes
:����������
random_color_affine_2/AddAddV2inputs(random_color_affine_2/random_uniform:z:0*
T0*/
_output_shapes
:���������
'p
.random_color_affine_2/random_uniform_1/shape/1Const*
_output_shapes
: *
dtype0*
value	B :p
.random_color_affine_2/random_uniform_1/shape/2Const*
_output_shapes
: *
dtype0*
value	B :p
.random_color_affine_2/random_uniform_1/shape/3Const*
_output_shapes
: *
dtype0*
value	B :�
,random_color_affine_2/random_uniform_1/shapePack,random_color_affine_2/strided_slice:output:07random_color_affine_2/random_uniform_1/shape/1:output:07random_color_affine_2/random_uniform_1/shape/2:output:07random_color_affine_2/random_uniform_1/shape/3:output:0*
N*
T0*
_output_shapes
:o
*random_color_affine_2/random_uniform_1/minConst*
_output_shapes
: *
dtype0*
valueB
 *��̽o
*random_color_affine_2/random_uniform_1/maxConst*
_output_shapes
: *
dtype0*
valueB
 *���=�
4random_color_affine_2/random_uniform_1/RandomUniformRandomUniform5random_color_affine_2/random_uniform_1/shape:output:0*
T0*/
_output_shapes
:���������*
dtype0�
*random_color_affine_2/random_uniform_1/subSub3random_color_affine_2/random_uniform_1/max:output:03random_color_affine_2/random_uniform_1/min:output:0*
T0*
_output_shapes
: �
*random_color_affine_2/random_uniform_1/mulMul=random_color_affine_2/random_uniform_1/RandomUniform:output:0.random_color_affine_2/random_uniform_1/sub:z:0*
T0*/
_output_shapes
:����������
&random_color_affine_2/random_uniform_1AddV2.random_color_affine_2/random_uniform_1/mul:z:03random_color_affine_2/random_uniform_1/min:output:0*
T0*/
_output_shapes
:����������
,random_color_affine_2/Mean/reduction_indicesConst*
_output_shapes
:*
dtype0*!
valueB"         �
random_color_affine_2/MeanMeanrandom_color_affine_2/Add:z:05random_color_affine_2/Mean/reduction_indices:output:0*
T0*/
_output_shapes
:���������*
	keep_dims(�
random_color_affine_2/subSubrandom_color_affine_2/Add:z:0#random_color_affine_2/Mean:output:0*
T0*/
_output_shapes
:���������
'�
random_color_affine_2/mulMul*random_color_affine_2/random_uniform_1:z:0random_color_affine_2/sub:z:0*
T0*/
_output_shapes
:���������
'�
random_color_affine_2/Add_1AddV2random_color_affine_2/Add:z:0random_color_affine_2/mul:z:0*
T0*/
_output_shapes
:���������
'r
-random_color_affine_2/clip_by_value/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?�
+random_color_affine_2/clip_by_value/MinimumMinimumrandom_color_affine_2/Add_1:z:06random_color_affine_2/clip_by_value/Minimum/y:output:0*
T0*/
_output_shapes
:���������
'j
%random_color_affine_2/clip_by_value/yConst*
_output_shapes
: *
dtype0*
valueB
 *    �
#random_color_affine_2/clip_by_valueMaximum/random_color_affine_2/clip_by_value/Minimum:z:0.random_color_affine_2/clip_by_value/y:output:0*
T0*/
_output_shapes
:���������
'w
IdentityIdentity'random_color_affine_2/clip_by_value:z:0*
T0*/
_output_shapes
:���������
'"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
':W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�
g
.__inference_gaussian_noise_layer_call_fn_11688

inputs
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_gaussian_noise_layer_call_and_return_conditional_losses_10460w
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:���������
'`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
'22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�#
o
P__inference_random_color_affine_2_layer_call_and_return_conditional_losses_10391

images
identity�I
ShapeShapeimages*
T0*
_output_shapes
::��]
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: _
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:_
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:�
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_maskX
random_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :X
random_uniform/shape/2Const*
_output_shapes
: *
dtype0*
value	B :X
random_uniform/shape/3Const*
_output_shapes
: *
dtype0*
value	B :�
random_uniform/shapePackstrided_slice:output:0random_uniform/shape/1:output:0random_uniform/shape/2:output:0random_uniform/shape/3:output:0*
N*
T0*
_output_shapes
:W
random_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *����W
random_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *���>�
random_uniform/RandomUniformRandomUniformrandom_uniform/shape:output:0*
T0*/
_output_shapes
:���������*
dtype0t
random_uniform/subSubrandom_uniform/max:output:0random_uniform/min:output:0*
T0*
_output_shapes
: �
random_uniform/mulMul%random_uniform/RandomUniform:output:0random_uniform/sub:z:0*
T0*/
_output_shapes
:����������
random_uniformAddV2random_uniform/mul:z:0random_uniform/min:output:0*
T0*/
_output_shapes
:���������b
AddAddV2imagesrandom_uniform:z:0*
T0*/
_output_shapes
:���������
'Z
random_uniform_1/shape/1Const*
_output_shapes
: *
dtype0*
value	B :Z
random_uniform_1/shape/2Const*
_output_shapes
: *
dtype0*
value	B :Z
random_uniform_1/shape/3Const*
_output_shapes
: *
dtype0*
value	B :�
random_uniform_1/shapePackstrided_slice:output:0!random_uniform_1/shape/1:output:0!random_uniform_1/shape/2:output:0!random_uniform_1/shape/3:output:0*
N*
T0*
_output_shapes
:Y
random_uniform_1/minConst*
_output_shapes
: *
dtype0*
valueB
 *��̽Y
random_uniform_1/maxConst*
_output_shapes
: *
dtype0*
valueB
 *���=�
random_uniform_1/RandomUniformRandomUniformrandom_uniform_1/shape:output:0*
T0*/
_output_shapes
:���������*
dtype0z
random_uniform_1/subSubrandom_uniform_1/max:output:0random_uniform_1/min:output:0*
T0*
_output_shapes
: �
random_uniform_1/mulMul'random_uniform_1/RandomUniform:output:0random_uniform_1/sub:z:0*
T0*/
_output_shapes
:����������
random_uniform_1AddV2random_uniform_1/mul:z:0random_uniform_1/min:output:0*
T0*/
_output_shapes
:���������k
Mean/reduction_indicesConst*
_output_shapes
:*
dtype0*!
valueB"         �
MeanMeanAdd:z:0Mean/reduction_indices:output:0*
T0*/
_output_shapes
:���������*
	keep_dims(\
subSubAdd:z:0Mean:output:0*
T0*/
_output_shapes
:���������
'c
mulMulrandom_uniform_1:z:0sub:z:0*
T0*/
_output_shapes
:���������
'Z
Add_1AddV2Add:z:0mul:z:0*
T0*/
_output_shapes
:���������
'\
clip_by_value/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?�
clip_by_value/MinimumMinimum	Add_1:z:0 clip_by_value/Minimum/y:output:0*
T0*/
_output_shapes
:���������
'T
clip_by_value/yConst*
_output_shapes
: *
dtype0*
valueB
 *    �
clip_by_valueMaximumclip_by_value/Minimum:z:0clip_by_value/y:output:0*
T0*/
_output_shapes
:���������
'a
IdentityIdentityclip_by_value:z:0*
T0*/
_output_shapes
:���������
'"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
':W S
/
_output_shapes
:���������
'
 
_user_specified_nameimages
Ս
�	
K__inference_finetuning_model_layer_call_and_return_conditional_losses_11334

inputsL
2sequential_2_conv_0_conv2d_readvariableop_resource: A
3sequential_2_conv_0_biasadd_readvariableop_resource: L
2sequential_2_conv_1_conv2d_readvariableop_resource:  A
3sequential_2_conv_1_biasadd_readvariableop_resource: C
0sequential_2_fc_0_matmul_readvariableop_resource:	�. ?
1sequential_2_fc_0_biasadd_readvariableop_resource: 8
&dense_3_matmul_readvariableop_resource: 
5
'dense_3_biasadd_readvariableop_resource:

identity��-conv_0/bias/Regularizer/L2Loss/ReadVariableOp�/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp�-conv_1/bias/Regularizer/L2Loss/ReadVariableOp�/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp�dense_3/BiasAdd/ReadVariableOp�dense_3/MatMul/ReadVariableOp�+fc_0/bias/Regularizer/L2Loss/ReadVariableOp�-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp�*sequential_2/conv_0/BiasAdd/ReadVariableOp�)sequential_2/conv_0/Conv2D/ReadVariableOp�*sequential_2/conv_1/BiasAdd/ReadVariableOp�)sequential_2/conv_1/Conv2D/ReadVariableOp�(sequential_2/fc_0/BiasAdd/ReadVariableOp�'sequential_2/fc_0/MatMul/ReadVariableOpl
(sequential_3/random_color_affine_2/ShapeShapeinputs*
T0*
_output_shapes
::���
6sequential_3/random_color_affine_2/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: �
8sequential_3/random_color_affine_2/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:�
8sequential_3/random_color_affine_2/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:�
0sequential_3/random_color_affine_2/strided_sliceStridedSlice1sequential_3/random_color_affine_2/Shape:output:0?sequential_3/random_color_affine_2/strided_slice/stack:output:0Asequential_3/random_color_affine_2/strided_slice/stack_1:output:0Asequential_3/random_color_affine_2/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask{
9sequential_3/random_color_affine_2/random_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :{
9sequential_3/random_color_affine_2/random_uniform/shape/2Const*
_output_shapes
: *
dtype0*
value	B :{
9sequential_3/random_color_affine_2/random_uniform/shape/3Const*
_output_shapes
: *
dtype0*
value	B :�
7sequential_3/random_color_affine_2/random_uniform/shapePack9sequential_3/random_color_affine_2/strided_slice:output:0Bsequential_3/random_color_affine_2/random_uniform/shape/1:output:0Bsequential_3/random_color_affine_2/random_uniform/shape/2:output:0Bsequential_3/random_color_affine_2/random_uniform/shape/3:output:0*
N*
T0*
_output_shapes
:z
5sequential_3/random_color_affine_2/random_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *����z
5sequential_3/random_color_affine_2/random_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *���>�
?sequential_3/random_color_affine_2/random_uniform/RandomUniformRandomUniform@sequential_3/random_color_affine_2/random_uniform/shape:output:0*
T0*/
_output_shapes
:���������*
dtype0�
5sequential_3/random_color_affine_2/random_uniform/subSub>sequential_3/random_color_affine_2/random_uniform/max:output:0>sequential_3/random_color_affine_2/random_uniform/min:output:0*
T0*
_output_shapes
: �
5sequential_3/random_color_affine_2/random_uniform/mulMulHsequential_3/random_color_affine_2/random_uniform/RandomUniform:output:09sequential_3/random_color_affine_2/random_uniform/sub:z:0*
T0*/
_output_shapes
:����������
1sequential_3/random_color_affine_2/random_uniformAddV29sequential_3/random_color_affine_2/random_uniform/mul:z:0>sequential_3/random_color_affine_2/random_uniform/min:output:0*
T0*/
_output_shapes
:����������
&sequential_3/random_color_affine_2/AddAddV2inputs5sequential_3/random_color_affine_2/random_uniform:z:0*
T0*/
_output_shapes
:���������
'}
;sequential_3/random_color_affine_2/random_uniform_1/shape/1Const*
_output_shapes
: *
dtype0*
value	B :}
;sequential_3/random_color_affine_2/random_uniform_1/shape/2Const*
_output_shapes
: *
dtype0*
value	B :}
;sequential_3/random_color_affine_2/random_uniform_1/shape/3Const*
_output_shapes
: *
dtype0*
value	B :�
9sequential_3/random_color_affine_2/random_uniform_1/shapePack9sequential_3/random_color_affine_2/strided_slice:output:0Dsequential_3/random_color_affine_2/random_uniform_1/shape/1:output:0Dsequential_3/random_color_affine_2/random_uniform_1/shape/2:output:0Dsequential_3/random_color_affine_2/random_uniform_1/shape/3:output:0*
N*
T0*
_output_shapes
:|
7sequential_3/random_color_affine_2/random_uniform_1/minConst*
_output_shapes
: *
dtype0*
valueB
 *��̽|
7sequential_3/random_color_affine_2/random_uniform_1/maxConst*
_output_shapes
: *
dtype0*
valueB
 *���=�
Asequential_3/random_color_affine_2/random_uniform_1/RandomUniformRandomUniformBsequential_3/random_color_affine_2/random_uniform_1/shape:output:0*
T0*/
_output_shapes
:���������*
dtype0�
7sequential_3/random_color_affine_2/random_uniform_1/subSub@sequential_3/random_color_affine_2/random_uniform_1/max:output:0@sequential_3/random_color_affine_2/random_uniform_1/min:output:0*
T0*
_output_shapes
: �
7sequential_3/random_color_affine_2/random_uniform_1/mulMulJsequential_3/random_color_affine_2/random_uniform_1/RandomUniform:output:0;sequential_3/random_color_affine_2/random_uniform_1/sub:z:0*
T0*/
_output_shapes
:����������
3sequential_3/random_color_affine_2/random_uniform_1AddV2;sequential_3/random_color_affine_2/random_uniform_1/mul:z:0@sequential_3/random_color_affine_2/random_uniform_1/min:output:0*
T0*/
_output_shapes
:����������
9sequential_3/random_color_affine_2/Mean/reduction_indicesConst*
_output_shapes
:*
dtype0*!
valueB"         �
'sequential_3/random_color_affine_2/MeanMean*sequential_3/random_color_affine_2/Add:z:0Bsequential_3/random_color_affine_2/Mean/reduction_indices:output:0*
T0*/
_output_shapes
:���������*
	keep_dims(�
&sequential_3/random_color_affine_2/subSub*sequential_3/random_color_affine_2/Add:z:00sequential_3/random_color_affine_2/Mean:output:0*
T0*/
_output_shapes
:���������
'�
&sequential_3/random_color_affine_2/mulMul7sequential_3/random_color_affine_2/random_uniform_1:z:0*sequential_3/random_color_affine_2/sub:z:0*
T0*/
_output_shapes
:���������
'�
(sequential_3/random_color_affine_2/Add_1AddV2*sequential_3/random_color_affine_2/Add:z:0*sequential_3/random_color_affine_2/mul:z:0*
T0*/
_output_shapes
:���������
'
:sequential_3/random_color_affine_2/clip_by_value/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?�
8sequential_3/random_color_affine_2/clip_by_value/MinimumMinimum,sequential_3/random_color_affine_2/Add_1:z:0Csequential_3/random_color_affine_2/clip_by_value/Minimum/y:output:0*
T0*/
_output_shapes
:���������
'w
2sequential_3/random_color_affine_2/clip_by_value/yConst*
_output_shapes
: *
dtype0*
valueB
 *    �
0sequential_3/random_color_affine_2/clip_by_valueMaximum<sequential_3/random_color_affine_2/clip_by_value/Minimum:z:0;sequential_3/random_color_affine_2/clip_by_value/y:output:0*
T0*/
_output_shapes
:���������
'�
!sequential_2/gaussian_noise/ShapeShape4sequential_3/random_color_affine_2/clip_by_value:z:0*
T0*
_output_shapes
::��s
.sequential_2/gaussian_noise/random_normal/meanConst*
_output_shapes
: *
dtype0*
valueB
 *    u
0sequential_2/gaussian_noise/random_normal/stddevConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
>sequential_2/gaussian_noise/random_normal/RandomStandardNormalRandomStandardNormal*sequential_2/gaussian_noise/Shape:output:0*
T0*/
_output_shapes
:���������
'*
dtype0�
-sequential_2/gaussian_noise/random_normal/mulMulGsequential_2/gaussian_noise/random_normal/RandomStandardNormal:output:09sequential_2/gaussian_noise/random_normal/stddev:output:0*
T0*/
_output_shapes
:���������
'�
)sequential_2/gaussian_noise/random_normalAddV21sequential_2/gaussian_noise/random_normal/mul:z:07sequential_2/gaussian_noise/random_normal/mean:output:0*
T0*/
_output_shapes
:���������
'�
sequential_2/gaussian_noise/addAddV24sequential_3/random_color_affine_2/clip_by_value:z:0-sequential_2/gaussian_noise/random_normal:z:0*
T0*/
_output_shapes
:���������
'�
)sequential_2/conv_0/Conv2D/ReadVariableOpReadVariableOp2sequential_2_conv_0_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype0�
sequential_2/conv_0/Conv2DConv2D#sequential_2/gaussian_noise/add:z:01sequential_2/conv_0/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������# *
paddingVALID*
strides
�
*sequential_2/conv_0/BiasAdd/ReadVariableOpReadVariableOp3sequential_2_conv_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
sequential_2/conv_0/BiasAddBiasAdd#sequential_2/conv_0/Conv2D:output:02sequential_2/conv_0/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������# �
sequential_2/conv_0/ReluRelu$sequential_2/conv_0/BiasAdd:output:0*
T0*/
_output_shapes
:���������# �
)sequential_2/conv_1/Conv2D/ReadVariableOpReadVariableOp2sequential_2_conv_1_conv2d_readvariableop_resource*&
_output_shapes
:  *
dtype0�
sequential_2/conv_1/Conv2DConv2D&sequential_2/conv_0/Relu:activations:01sequential_2/conv_1/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:��������� *
paddingVALID*
strides
�
*sequential_2/conv_1/BiasAdd/ReadVariableOpReadVariableOp3sequential_2_conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
sequential_2/conv_1/BiasAddBiasAdd#sequential_2/conv_1/Conv2D:output:02sequential_2/conv_1/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:��������� �
sequential_2/conv_1/ReluRelu$sequential_2/conv_1/BiasAdd:output:0*
T0*/
_output_shapes
:��������� k
sequential_2/flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"����@  �
sequential_2/flatten/ReshapeReshape&sequential_2/conv_1/Relu:activations:0#sequential_2/flatten/Const:output:0*
T0*(
_output_shapes
:����������.�
'sequential_2/fc_0/MatMul/ReadVariableOpReadVariableOp0sequential_2_fc_0_matmul_readvariableop_resource*
_output_shapes
:	�. *
dtype0�
sequential_2/fc_0/MatMulMatMul%sequential_2/flatten/Reshape:output:0/sequential_2/fc_0/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� �
(sequential_2/fc_0/BiasAdd/ReadVariableOpReadVariableOp1sequential_2_fc_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
sequential_2/fc_0/BiasAddBiasAdd"sequential_2/fc_0/MatMul:product:00sequential_2/fc_0/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� t
sequential_2/fc_0/ReluRelu"sequential_2/fc_0/BiasAdd:output:0*
T0*'
_output_shapes
:��������� �
dense_3/MatMul/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource*
_output_shapes

: 
*
dtype0�
dense_3/MatMulMatMul$sequential_2/fc_0/Relu:activations:0%dense_3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
dense_3/BiasAdd/ReadVariableOpReadVariableOp'dense_3_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype0�
dense_3/BiasAddBiasAdddense_3/MatMul:product:0&dense_3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp2sequential_2_conv_0_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype0�
 conv_0/kernel/Regularizer/L2LossL2Loss7conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/kernel/Regularizer/mulMul(conv_0/kernel/Regularizer/mul/x:output:0)conv_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-conv_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp3sequential_2_conv_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_0/bias/Regularizer/L2LossL2Loss5conv_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/bias/Regularizer/mulMul&conv_0/bias/Regularizer/mul/x:output:0'conv_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp2sequential_2_conv_1_conv2d_readvariableop_resource*&
_output_shapes
:  *
dtype0�
 conv_1/kernel/Regularizer/L2LossL2Loss7conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/kernel/Regularizer/mulMul(conv_1/kernel/Regularizer/mul/x:output:0)conv_1/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-conv_1/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp3sequential_2_conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_1/bias/Regularizer/L2LossL2Loss5conv_1/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/bias/Regularizer/mulMul&conv_1/bias/Regularizer/mul/x:output:0'conv_1/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp0sequential_2_fc_0_matmul_readvariableop_resource*
_output_shapes
:	�. *
dtype0�
fc_0/kernel/Regularizer/L2LossL2Loss5fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
fc_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/kernel/Regularizer/mulMul&fc_0/kernel/Regularizer/mul/x:output:0'fc_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
+fc_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp1sequential_2_fc_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0|
fc_0/bias/Regularizer/L2LossL2Loss3fc_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: `
fc_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/bias/Regularizer/mulMul$fc_0/bias/Regularizer/mul/x:output:0%fc_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: g
IdentityIdentitydense_3/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������
�
NoOpNoOp.^conv_0/bias/Regularizer/L2Loss/ReadVariableOp0^conv_0/kernel/Regularizer/L2Loss/ReadVariableOp.^conv_1/bias/Regularizer/L2Loss/ReadVariableOp0^conv_1/kernel/Regularizer/L2Loss/ReadVariableOp^dense_3/BiasAdd/ReadVariableOp^dense_3/MatMul/ReadVariableOp,^fc_0/bias/Regularizer/L2Loss/ReadVariableOp.^fc_0/kernel/Regularizer/L2Loss/ReadVariableOp+^sequential_2/conv_0/BiasAdd/ReadVariableOp*^sequential_2/conv_0/Conv2D/ReadVariableOp+^sequential_2/conv_1/BiasAdd/ReadVariableOp*^sequential_2/conv_1/Conv2D/ReadVariableOp)^sequential_2/fc_0/BiasAdd/ReadVariableOp(^sequential_2/fc_0/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:���������
': : : : : : : : 2^
-conv_0/bias/Regularizer/L2Loss/ReadVariableOp-conv_0/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp2^
-conv_1/bias/Regularizer/L2Loss/ReadVariableOp-conv_1/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp2@
dense_3/BiasAdd/ReadVariableOpdense_3/BiasAdd/ReadVariableOp2>
dense_3/MatMul/ReadVariableOpdense_3/MatMul/ReadVariableOp2Z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOp+fc_0/bias/Regularizer/L2Loss/ReadVariableOp2^
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp2X
*sequential_2/conv_0/BiasAdd/ReadVariableOp*sequential_2/conv_0/BiasAdd/ReadVariableOp2V
)sequential_2/conv_0/Conv2D/ReadVariableOp)sequential_2/conv_0/Conv2D/ReadVariableOp2X
*sequential_2/conv_1/BiasAdd/ReadVariableOp*sequential_2/conv_1/BiasAdd/ReadVariableOp2V
)sequential_2/conv_1/Conv2D/ReadVariableOp)sequential_2/conv_1/Conv2D/ReadVariableOp2T
(sequential_2/fc_0/BiasAdd/ReadVariableOp(sequential_2/fc_0/BiasAdd/ReadVariableOp2R
'sequential_2/fc_0/MatMul/ReadVariableOp'sequential_2/fc_0/MatMul/ReadVariableOp:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�	
�
__inference_loss_fn_0_11812R
8conv_0_kernel_regularizer_l2loss_readvariableop_resource: 
identity��/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp�
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp8conv_0_kernel_regularizer_l2loss_readvariableop_resource*&
_output_shapes
: *
dtype0�
 conv_0/kernel/Regularizer/L2LossL2Loss7conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/kernel/Regularizer/mulMul(conv_0/kernel/Regularizer/mul/x:output:0)conv_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: _
IdentityIdentity!conv_0/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: x
NoOpNoOp0^conv_0/kernel/Regularizer/L2Loss/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2b
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp
�	
�
#__inference_signature_wrapper_11169
input_5!
unknown: 
	unknown_0: #
	unknown_1:  
	unknown_2: 
	unknown_3:	�. 
	unknown_4: 
	unknown_5: 

	unknown_6:

identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_5unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *)
f$R"
 __inference__wrapped_model_10348o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������
`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:���������
': : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:X T
/
_output_shapes
:���������
'
!
_user_specified_name	input_5
�5
�
K__inference_finetuning_model_layer_call_and_return_conditional_losses_10981

inputs,
sequential_2_10938:  
sequential_2_10940: ,
sequential_2_10942:   
sequential_2_10944: %
sequential_2_10946:	�.  
sequential_2_10948: 
dense_3_10951: 

dense_3_10953:

identity��-conv_0/bias/Regularizer/L2Loss/ReadVariableOp�/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp�-conv_1/bias/Regularizer/L2Loss/ReadVariableOp�/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp�dense_3/StatefulPartitionedCall�+fc_0/bias/Regularizer/L2Loss/ReadVariableOp�-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp�$sequential_2/StatefulPartitionedCall�$sequential_3/StatefulPartitionedCall�
$sequential_3/StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_3_layer_call_and_return_conditional_losses_10411�
$sequential_2/StatefulPartitionedCallStatefulPartitionedCall-sequential_3/StatefulPartitionedCall:output:0sequential_2_10938sequential_2_10940sequential_2_10942sequential_2_10944sequential_2_10946sequential_2_10948*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *(
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_2_layer_call_and_return_conditional_losses_10667�
dense_3/StatefulPartitionedCallStatefulPartitionedCall-sequential_2/StatefulPartitionedCall:output:0dense_3_10951dense_3_10953*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_3_layer_call_and_return_conditional_losses_10853�
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10938*&
_output_shapes
: *
dtype0�
 conv_0/kernel/Regularizer/L2LossL2Loss7conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/kernel/Regularizer/mulMul(conv_0/kernel/Regularizer/mul/x:output:0)conv_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: |
-conv_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10940*
_output_shapes
: *
dtype0�
conv_0/bias/Regularizer/L2LossL2Loss5conv_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/bias/Regularizer/mulMul&conv_0/bias/Regularizer/mul/x:output:0'conv_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10942*&
_output_shapes
:  *
dtype0�
 conv_1/kernel/Regularizer/L2LossL2Loss7conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/kernel/Regularizer/mulMul(conv_1/kernel/Regularizer/mul/x:output:0)conv_1/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: |
-conv_1/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10944*
_output_shapes
: *
dtype0�
conv_1/bias/Regularizer/L2LossL2Loss5conv_1/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/bias/Regularizer/mulMul&conv_1/bias/Regularizer/mul/x:output:0'conv_1/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10946*
_output_shapes
:	�. *
dtype0�
fc_0/kernel/Regularizer/L2LossL2Loss5fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
fc_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/kernel/Regularizer/mulMul&fc_0/kernel/Regularizer/mul/x:output:0'fc_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10948*
_output_shapes
: *
dtype0|
fc_0/bias/Regularizer/L2LossL2Loss3fc_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: `
fc_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/bias/Regularizer/mulMul$fc_0/bias/Regularizer/mul/x:output:0%fc_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: w
IdentityIdentity(dense_3/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������
�
NoOpNoOp.^conv_0/bias/Regularizer/L2Loss/ReadVariableOp0^conv_0/kernel/Regularizer/L2Loss/ReadVariableOp.^conv_1/bias/Regularizer/L2Loss/ReadVariableOp0^conv_1/kernel/Regularizer/L2Loss/ReadVariableOp ^dense_3/StatefulPartitionedCall,^fc_0/bias/Regularizer/L2Loss/ReadVariableOp.^fc_0/kernel/Regularizer/L2Loss/ReadVariableOp%^sequential_2/StatefulPartitionedCall%^sequential_3/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:���������
': : : : : : : : 2^
-conv_0/bias/Regularizer/L2Loss/ReadVariableOp-conv_0/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp2^
-conv_1/bias/Regularizer/L2Loss/ReadVariableOp-conv_1/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2Z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOp+fc_0/bias/Regularizer/L2Loss/ReadVariableOp2^
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp2L
$sequential_2/StatefulPartitionedCall$sequential_2/StatefulPartitionedCall2L
$sequential_3/StatefulPartitionedCall$sequential_3/StatefulPartitionedCall:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�
�
__inference_loss_fn_3_11839D
6conv_1_bias_regularizer_l2loss_readvariableop_resource: 
identity��-conv_1/bias/Regularizer/L2Loss/ReadVariableOp�
-conv_1/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp6conv_1_bias_regularizer_l2loss_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_1/bias/Regularizer/L2LossL2Loss5conv_1/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/bias/Regularizer/mulMul&conv_1/bias/Regularizer/mul/x:output:0'conv_1/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: ]
IdentityIdentityconv_1/bias/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: v
NoOpNoOp.^conv_1/bias/Regularizer/L2Loss/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2^
-conv_1/bias/Regularizer/L2Loss/ReadVariableOp-conv_1/bias/Regularizer/L2Loss/ReadVariableOp
�
Q
5__inference_random_color_affine_2_layer_call_fn_11640

images
identity�
PartitionedCallPartitionedCallimages*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *Y
fTRR
P__inference_random_color_affine_2_layer_call_and_return_conditional_losses_10400h
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:���������
'"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
':W S
/
_output_shapes
:���������
'
 
_user_specified_nameimages
��
�
__inference__traced_save_12060
file_prefix7
%read_disablecopyonread_dense_3_kernel: 
3
%read_1_disablecopyonread_dense_3_bias:
@
&read_2_disablecopyonread_conv_0_kernel: 2
$read_3_disablecopyonread_conv_0_bias: @
&read_4_disablecopyonread_conv_1_kernel:  2
$read_5_disablecopyonread_conv_1_bias: 7
$read_6_disablecopyonread_fc_0_kernel:	�. 0
"read_7_disablecopyonread_fc_0_bias: ,
"read_8_disablecopyonread_iteration:	 0
&read_9_disablecopyonread_learning_rate: H
.read_10_disablecopyonread_adam_m_conv_0_kernel: H
.read_11_disablecopyonread_adam_v_conv_0_kernel: :
,read_12_disablecopyonread_adam_m_conv_0_bias: :
,read_13_disablecopyonread_adam_v_conv_0_bias: H
.read_14_disablecopyonread_adam_m_conv_1_kernel:  H
.read_15_disablecopyonread_adam_v_conv_1_kernel:  :
,read_16_disablecopyonread_adam_m_conv_1_bias: :
,read_17_disablecopyonread_adam_v_conv_1_bias: ?
,read_18_disablecopyonread_adam_m_fc_0_kernel:	�. ?
,read_19_disablecopyonread_adam_v_fc_0_kernel:	�. 8
*read_20_disablecopyonread_adam_m_fc_0_bias: 8
*read_21_disablecopyonread_adam_v_fc_0_bias: A
/read_22_disablecopyonread_adam_m_dense_3_kernel: 
A
/read_23_disablecopyonread_adam_v_dense_3_kernel: 
;
-read_24_disablecopyonread_adam_m_dense_3_bias:
;
-read_25_disablecopyonread_adam_v_dense_3_bias:
+
!read_26_disablecopyonread_total_1: +
!read_27_disablecopyonread_count_1: )
read_28_disablecopyonread_total: )
read_29_disablecopyonread_count: 
savev2_const
identity_61��MergeV2Checkpoints�Read/DisableCopyOnRead�Read/ReadVariableOp�Read_1/DisableCopyOnRead�Read_1/ReadVariableOp�Read_10/DisableCopyOnRead�Read_10/ReadVariableOp�Read_11/DisableCopyOnRead�Read_11/ReadVariableOp�Read_12/DisableCopyOnRead�Read_12/ReadVariableOp�Read_13/DisableCopyOnRead�Read_13/ReadVariableOp�Read_14/DisableCopyOnRead�Read_14/ReadVariableOp�Read_15/DisableCopyOnRead�Read_15/ReadVariableOp�Read_16/DisableCopyOnRead�Read_16/ReadVariableOp�Read_17/DisableCopyOnRead�Read_17/ReadVariableOp�Read_18/DisableCopyOnRead�Read_18/ReadVariableOp�Read_19/DisableCopyOnRead�Read_19/ReadVariableOp�Read_2/DisableCopyOnRead�Read_2/ReadVariableOp�Read_20/DisableCopyOnRead�Read_20/ReadVariableOp�Read_21/DisableCopyOnRead�Read_21/ReadVariableOp�Read_22/DisableCopyOnRead�Read_22/ReadVariableOp�Read_23/DisableCopyOnRead�Read_23/ReadVariableOp�Read_24/DisableCopyOnRead�Read_24/ReadVariableOp�Read_25/DisableCopyOnRead�Read_25/ReadVariableOp�Read_26/DisableCopyOnRead�Read_26/ReadVariableOp�Read_27/DisableCopyOnRead�Read_27/ReadVariableOp�Read_28/DisableCopyOnRead�Read_28/ReadVariableOp�Read_29/DisableCopyOnRead�Read_29/ReadVariableOp�Read_3/DisableCopyOnRead�Read_3/ReadVariableOp�Read_4/DisableCopyOnRead�Read_4/ReadVariableOp�Read_5/DisableCopyOnRead�Read_5/ReadVariableOp�Read_6/DisableCopyOnRead�Read_6/ReadVariableOp�Read_7/DisableCopyOnRead�Read_7/ReadVariableOp�Read_8/DisableCopyOnRead�Read_8/ReadVariableOp�Read_9/DisableCopyOnRead�Read_9/ReadVariableOpw
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*Z
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.parta
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B
_temp/part�
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: f

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: L

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :f
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : �
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: w
Read/DisableCopyOnReadDisableCopyOnRead%read_disablecopyonread_dense_3_kernel"/device:CPU:0*
_output_shapes
 �
Read/ReadVariableOpReadVariableOp%read_disablecopyonread_dense_3_kernel^Read/DisableCopyOnRead"/device:CPU:0*
_output_shapes

: 
*
dtype0i
IdentityIdentityRead/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

: 
a

Identity_1IdentityIdentity:output:0"/device:CPU:0*
T0*
_output_shapes

: 
y
Read_1/DisableCopyOnReadDisableCopyOnRead%read_1_disablecopyonread_dense_3_bias"/device:CPU:0*
_output_shapes
 �
Read_1/ReadVariableOpReadVariableOp%read_1_disablecopyonread_dense_3_bias^Read_1/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:
*
dtype0i

Identity_2IdentityRead_1/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:
_

Identity_3IdentityIdentity_2:output:0"/device:CPU:0*
T0*
_output_shapes
:
z
Read_2/DisableCopyOnReadDisableCopyOnRead&read_2_disablecopyonread_conv_0_kernel"/device:CPU:0*
_output_shapes
 �
Read_2/ReadVariableOpReadVariableOp&read_2_disablecopyonread_conv_0_kernel^Read_2/DisableCopyOnRead"/device:CPU:0*&
_output_shapes
: *
dtype0u

Identity_4IdentityRead_2/ReadVariableOp:value:0"/device:CPU:0*
T0*&
_output_shapes
: k

Identity_5IdentityIdentity_4:output:0"/device:CPU:0*
T0*&
_output_shapes
: x
Read_3/DisableCopyOnReadDisableCopyOnRead$read_3_disablecopyonread_conv_0_bias"/device:CPU:0*
_output_shapes
 �
Read_3/ReadVariableOpReadVariableOp$read_3_disablecopyonread_conv_0_bias^Read_3/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0i

Identity_6IdentityRead_3/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: _

Identity_7IdentityIdentity_6:output:0"/device:CPU:0*
T0*
_output_shapes
: z
Read_4/DisableCopyOnReadDisableCopyOnRead&read_4_disablecopyonread_conv_1_kernel"/device:CPU:0*
_output_shapes
 �
Read_4/ReadVariableOpReadVariableOp&read_4_disablecopyonread_conv_1_kernel^Read_4/DisableCopyOnRead"/device:CPU:0*&
_output_shapes
:  *
dtype0u

Identity_8IdentityRead_4/ReadVariableOp:value:0"/device:CPU:0*
T0*&
_output_shapes
:  k

Identity_9IdentityIdentity_8:output:0"/device:CPU:0*
T0*&
_output_shapes
:  x
Read_5/DisableCopyOnReadDisableCopyOnRead$read_5_disablecopyonread_conv_1_bias"/device:CPU:0*
_output_shapes
 �
Read_5/ReadVariableOpReadVariableOp$read_5_disablecopyonread_conv_1_bias^Read_5/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0j
Identity_10IdentityRead_5/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: a
Identity_11IdentityIdentity_10:output:0"/device:CPU:0*
T0*
_output_shapes
: x
Read_6/DisableCopyOnReadDisableCopyOnRead$read_6_disablecopyonread_fc_0_kernel"/device:CPU:0*
_output_shapes
 �
Read_6/ReadVariableOpReadVariableOp$read_6_disablecopyonread_fc_0_kernel^Read_6/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:	�. *
dtype0o
Identity_12IdentityRead_6/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:	�. f
Identity_13IdentityIdentity_12:output:0"/device:CPU:0*
T0*
_output_shapes
:	�. v
Read_7/DisableCopyOnReadDisableCopyOnRead"read_7_disablecopyonread_fc_0_bias"/device:CPU:0*
_output_shapes
 �
Read_7/ReadVariableOpReadVariableOp"read_7_disablecopyonread_fc_0_bias^Read_7/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0j
Identity_14IdentityRead_7/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: a
Identity_15IdentityIdentity_14:output:0"/device:CPU:0*
T0*
_output_shapes
: v
Read_8/DisableCopyOnReadDisableCopyOnRead"read_8_disablecopyonread_iteration"/device:CPU:0*
_output_shapes
 �
Read_8/ReadVariableOpReadVariableOp"read_8_disablecopyonread_iteration^Read_8/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0	f
Identity_16IdentityRead_8/ReadVariableOp:value:0"/device:CPU:0*
T0	*
_output_shapes
: ]
Identity_17IdentityIdentity_16:output:0"/device:CPU:0*
T0	*
_output_shapes
: z
Read_9/DisableCopyOnReadDisableCopyOnRead&read_9_disablecopyonread_learning_rate"/device:CPU:0*
_output_shapes
 �
Read_9/ReadVariableOpReadVariableOp&read_9_disablecopyonread_learning_rate^Read_9/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0f
Identity_18IdentityRead_9/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: ]
Identity_19IdentityIdentity_18:output:0"/device:CPU:0*
T0*
_output_shapes
: �
Read_10/DisableCopyOnReadDisableCopyOnRead.read_10_disablecopyonread_adam_m_conv_0_kernel"/device:CPU:0*
_output_shapes
 �
Read_10/ReadVariableOpReadVariableOp.read_10_disablecopyonread_adam_m_conv_0_kernel^Read_10/DisableCopyOnRead"/device:CPU:0*&
_output_shapes
: *
dtype0w
Identity_20IdentityRead_10/ReadVariableOp:value:0"/device:CPU:0*
T0*&
_output_shapes
: m
Identity_21IdentityIdentity_20:output:0"/device:CPU:0*
T0*&
_output_shapes
: �
Read_11/DisableCopyOnReadDisableCopyOnRead.read_11_disablecopyonread_adam_v_conv_0_kernel"/device:CPU:0*
_output_shapes
 �
Read_11/ReadVariableOpReadVariableOp.read_11_disablecopyonread_adam_v_conv_0_kernel^Read_11/DisableCopyOnRead"/device:CPU:0*&
_output_shapes
: *
dtype0w
Identity_22IdentityRead_11/ReadVariableOp:value:0"/device:CPU:0*
T0*&
_output_shapes
: m
Identity_23IdentityIdentity_22:output:0"/device:CPU:0*
T0*&
_output_shapes
: �
Read_12/DisableCopyOnReadDisableCopyOnRead,read_12_disablecopyonread_adam_m_conv_0_bias"/device:CPU:0*
_output_shapes
 �
Read_12/ReadVariableOpReadVariableOp,read_12_disablecopyonread_adam_m_conv_0_bias^Read_12/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0k
Identity_24IdentityRead_12/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: a
Identity_25IdentityIdentity_24:output:0"/device:CPU:0*
T0*
_output_shapes
: �
Read_13/DisableCopyOnReadDisableCopyOnRead,read_13_disablecopyonread_adam_v_conv_0_bias"/device:CPU:0*
_output_shapes
 �
Read_13/ReadVariableOpReadVariableOp,read_13_disablecopyonread_adam_v_conv_0_bias^Read_13/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0k
Identity_26IdentityRead_13/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: a
Identity_27IdentityIdentity_26:output:0"/device:CPU:0*
T0*
_output_shapes
: �
Read_14/DisableCopyOnReadDisableCopyOnRead.read_14_disablecopyonread_adam_m_conv_1_kernel"/device:CPU:0*
_output_shapes
 �
Read_14/ReadVariableOpReadVariableOp.read_14_disablecopyonread_adam_m_conv_1_kernel^Read_14/DisableCopyOnRead"/device:CPU:0*&
_output_shapes
:  *
dtype0w
Identity_28IdentityRead_14/ReadVariableOp:value:0"/device:CPU:0*
T0*&
_output_shapes
:  m
Identity_29IdentityIdentity_28:output:0"/device:CPU:0*
T0*&
_output_shapes
:  �
Read_15/DisableCopyOnReadDisableCopyOnRead.read_15_disablecopyonread_adam_v_conv_1_kernel"/device:CPU:0*
_output_shapes
 �
Read_15/ReadVariableOpReadVariableOp.read_15_disablecopyonread_adam_v_conv_1_kernel^Read_15/DisableCopyOnRead"/device:CPU:0*&
_output_shapes
:  *
dtype0w
Identity_30IdentityRead_15/ReadVariableOp:value:0"/device:CPU:0*
T0*&
_output_shapes
:  m
Identity_31IdentityIdentity_30:output:0"/device:CPU:0*
T0*&
_output_shapes
:  �
Read_16/DisableCopyOnReadDisableCopyOnRead,read_16_disablecopyonread_adam_m_conv_1_bias"/device:CPU:0*
_output_shapes
 �
Read_16/ReadVariableOpReadVariableOp,read_16_disablecopyonread_adam_m_conv_1_bias^Read_16/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0k
Identity_32IdentityRead_16/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: a
Identity_33IdentityIdentity_32:output:0"/device:CPU:0*
T0*
_output_shapes
: �
Read_17/DisableCopyOnReadDisableCopyOnRead,read_17_disablecopyonread_adam_v_conv_1_bias"/device:CPU:0*
_output_shapes
 �
Read_17/ReadVariableOpReadVariableOp,read_17_disablecopyonread_adam_v_conv_1_bias^Read_17/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0k
Identity_34IdentityRead_17/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: a
Identity_35IdentityIdentity_34:output:0"/device:CPU:0*
T0*
_output_shapes
: �
Read_18/DisableCopyOnReadDisableCopyOnRead,read_18_disablecopyonread_adam_m_fc_0_kernel"/device:CPU:0*
_output_shapes
 �
Read_18/ReadVariableOpReadVariableOp,read_18_disablecopyonread_adam_m_fc_0_kernel^Read_18/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:	�. *
dtype0p
Identity_36IdentityRead_18/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:	�. f
Identity_37IdentityIdentity_36:output:0"/device:CPU:0*
T0*
_output_shapes
:	�. �
Read_19/DisableCopyOnReadDisableCopyOnRead,read_19_disablecopyonread_adam_v_fc_0_kernel"/device:CPU:0*
_output_shapes
 �
Read_19/ReadVariableOpReadVariableOp,read_19_disablecopyonread_adam_v_fc_0_kernel^Read_19/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:	�. *
dtype0p
Identity_38IdentityRead_19/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:	�. f
Identity_39IdentityIdentity_38:output:0"/device:CPU:0*
T0*
_output_shapes
:	�. 
Read_20/DisableCopyOnReadDisableCopyOnRead*read_20_disablecopyonread_adam_m_fc_0_bias"/device:CPU:0*
_output_shapes
 �
Read_20/ReadVariableOpReadVariableOp*read_20_disablecopyonread_adam_m_fc_0_bias^Read_20/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0k
Identity_40IdentityRead_20/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: a
Identity_41IdentityIdentity_40:output:0"/device:CPU:0*
T0*
_output_shapes
: 
Read_21/DisableCopyOnReadDisableCopyOnRead*read_21_disablecopyonread_adam_v_fc_0_bias"/device:CPU:0*
_output_shapes
 �
Read_21/ReadVariableOpReadVariableOp*read_21_disablecopyonread_adam_v_fc_0_bias^Read_21/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0k
Identity_42IdentityRead_21/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: a
Identity_43IdentityIdentity_42:output:0"/device:CPU:0*
T0*
_output_shapes
: �
Read_22/DisableCopyOnReadDisableCopyOnRead/read_22_disablecopyonread_adam_m_dense_3_kernel"/device:CPU:0*
_output_shapes
 �
Read_22/ReadVariableOpReadVariableOp/read_22_disablecopyonread_adam_m_dense_3_kernel^Read_22/DisableCopyOnRead"/device:CPU:0*
_output_shapes

: 
*
dtype0o
Identity_44IdentityRead_22/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

: 
e
Identity_45IdentityIdentity_44:output:0"/device:CPU:0*
T0*
_output_shapes

: 
�
Read_23/DisableCopyOnReadDisableCopyOnRead/read_23_disablecopyonread_adam_v_dense_3_kernel"/device:CPU:0*
_output_shapes
 �
Read_23/ReadVariableOpReadVariableOp/read_23_disablecopyonread_adam_v_dense_3_kernel^Read_23/DisableCopyOnRead"/device:CPU:0*
_output_shapes

: 
*
dtype0o
Identity_46IdentityRead_23/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes

: 
e
Identity_47IdentityIdentity_46:output:0"/device:CPU:0*
T0*
_output_shapes

: 
�
Read_24/DisableCopyOnReadDisableCopyOnRead-read_24_disablecopyonread_adam_m_dense_3_bias"/device:CPU:0*
_output_shapes
 �
Read_24/ReadVariableOpReadVariableOp-read_24_disablecopyonread_adam_m_dense_3_bias^Read_24/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:
*
dtype0k
Identity_48IdentityRead_24/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:
a
Identity_49IdentityIdentity_48:output:0"/device:CPU:0*
T0*
_output_shapes
:
�
Read_25/DisableCopyOnReadDisableCopyOnRead-read_25_disablecopyonread_adam_v_dense_3_bias"/device:CPU:0*
_output_shapes
 �
Read_25/ReadVariableOpReadVariableOp-read_25_disablecopyonread_adam_v_dense_3_bias^Read_25/DisableCopyOnRead"/device:CPU:0*
_output_shapes
:
*
dtype0k
Identity_50IdentityRead_25/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
:
a
Identity_51IdentityIdentity_50:output:0"/device:CPU:0*
T0*
_output_shapes
:
v
Read_26/DisableCopyOnReadDisableCopyOnRead!read_26_disablecopyonread_total_1"/device:CPU:0*
_output_shapes
 �
Read_26/ReadVariableOpReadVariableOp!read_26_disablecopyonread_total_1^Read_26/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0g
Identity_52IdentityRead_26/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: ]
Identity_53IdentityIdentity_52:output:0"/device:CPU:0*
T0*
_output_shapes
: v
Read_27/DisableCopyOnReadDisableCopyOnRead!read_27_disablecopyonread_count_1"/device:CPU:0*
_output_shapes
 �
Read_27/ReadVariableOpReadVariableOp!read_27_disablecopyonread_count_1^Read_27/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0g
Identity_54IdentityRead_27/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: ]
Identity_55IdentityIdentity_54:output:0"/device:CPU:0*
T0*
_output_shapes
: t
Read_28/DisableCopyOnReadDisableCopyOnReadread_28_disablecopyonread_total"/device:CPU:0*
_output_shapes
 �
Read_28/ReadVariableOpReadVariableOpread_28_disablecopyonread_total^Read_28/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0g
Identity_56IdentityRead_28/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: ]
Identity_57IdentityIdentity_56:output:0"/device:CPU:0*
T0*
_output_shapes
: t
Read_29/DisableCopyOnReadDisableCopyOnReadread_29_disablecopyonread_count"/device:CPU:0*
_output_shapes
 �
Read_29/ReadVariableOpReadVariableOpread_29_disablecopyonread_count^Read_29/DisableCopyOnRead"/device:CPU:0*
_output_shapes
: *
dtype0g
Identity_58IdentityRead_29/ReadVariableOp:value:0"/device:CPU:0*
T0*
_output_shapes
: ]
Identity_59IdentityIdentity_58:output:0"/device:CPU:0*
T0*
_output_shapes
: �
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*�
value�B�B6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUEB3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*Q
valueHBFB B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B �
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0Identity_1:output:0Identity_3:output:0Identity_5:output:0Identity_7:output:0Identity_9:output:0Identity_11:output:0Identity_13:output:0Identity_15:output:0Identity_17:output:0Identity_19:output:0Identity_21:output:0Identity_23:output:0Identity_25:output:0Identity_27:output:0Identity_29:output:0Identity_31:output:0Identity_33:output:0Identity_35:output:0Identity_37:output:0Identity_39:output:0Identity_41:output:0Identity_43:output:0Identity_45:output:0Identity_47:output:0Identity_49:output:0Identity_51:output:0Identity_53:output:0Identity_55:output:0Identity_57:output:0Identity_59:output:0savev2_const"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *-
dtypes#
!2	�
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:�
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 i
Identity_60Identityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: U
Identity_61IdentityIdentity_60:output:0^NoOp*
T0*
_output_shapes
: �
NoOpNoOp^MergeV2Checkpoints^Read/DisableCopyOnRead^Read/ReadVariableOp^Read_1/DisableCopyOnRead^Read_1/ReadVariableOp^Read_10/DisableCopyOnRead^Read_10/ReadVariableOp^Read_11/DisableCopyOnRead^Read_11/ReadVariableOp^Read_12/DisableCopyOnRead^Read_12/ReadVariableOp^Read_13/DisableCopyOnRead^Read_13/ReadVariableOp^Read_14/DisableCopyOnRead^Read_14/ReadVariableOp^Read_15/DisableCopyOnRead^Read_15/ReadVariableOp^Read_16/DisableCopyOnRead^Read_16/ReadVariableOp^Read_17/DisableCopyOnRead^Read_17/ReadVariableOp^Read_18/DisableCopyOnRead^Read_18/ReadVariableOp^Read_19/DisableCopyOnRead^Read_19/ReadVariableOp^Read_2/DisableCopyOnRead^Read_2/ReadVariableOp^Read_20/DisableCopyOnRead^Read_20/ReadVariableOp^Read_21/DisableCopyOnRead^Read_21/ReadVariableOp^Read_22/DisableCopyOnRead^Read_22/ReadVariableOp^Read_23/DisableCopyOnRead^Read_23/ReadVariableOp^Read_24/DisableCopyOnRead^Read_24/ReadVariableOp^Read_25/DisableCopyOnRead^Read_25/ReadVariableOp^Read_26/DisableCopyOnRead^Read_26/ReadVariableOp^Read_27/DisableCopyOnRead^Read_27/ReadVariableOp^Read_28/DisableCopyOnRead^Read_28/ReadVariableOp^Read_29/DisableCopyOnRead^Read_29/ReadVariableOp^Read_3/DisableCopyOnRead^Read_3/ReadVariableOp^Read_4/DisableCopyOnRead^Read_4/ReadVariableOp^Read_5/DisableCopyOnRead^Read_5/ReadVariableOp^Read_6/DisableCopyOnRead^Read_6/ReadVariableOp^Read_7/DisableCopyOnRead^Read_7/ReadVariableOp^Read_8/DisableCopyOnRead^Read_8/ReadVariableOp^Read_9/DisableCopyOnRead^Read_9/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "#
identity_61Identity_61:output:0*S
_input_shapesB
@: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 2(
MergeV2CheckpointsMergeV2Checkpoints20
Read/DisableCopyOnReadRead/DisableCopyOnRead2*
Read/ReadVariableOpRead/ReadVariableOp24
Read_1/DisableCopyOnReadRead_1/DisableCopyOnRead2.
Read_1/ReadVariableOpRead_1/ReadVariableOp26
Read_10/DisableCopyOnReadRead_10/DisableCopyOnRead20
Read_10/ReadVariableOpRead_10/ReadVariableOp26
Read_11/DisableCopyOnReadRead_11/DisableCopyOnRead20
Read_11/ReadVariableOpRead_11/ReadVariableOp26
Read_12/DisableCopyOnReadRead_12/DisableCopyOnRead20
Read_12/ReadVariableOpRead_12/ReadVariableOp26
Read_13/DisableCopyOnReadRead_13/DisableCopyOnRead20
Read_13/ReadVariableOpRead_13/ReadVariableOp26
Read_14/DisableCopyOnReadRead_14/DisableCopyOnRead20
Read_14/ReadVariableOpRead_14/ReadVariableOp26
Read_15/DisableCopyOnReadRead_15/DisableCopyOnRead20
Read_15/ReadVariableOpRead_15/ReadVariableOp26
Read_16/DisableCopyOnReadRead_16/DisableCopyOnRead20
Read_16/ReadVariableOpRead_16/ReadVariableOp26
Read_17/DisableCopyOnReadRead_17/DisableCopyOnRead20
Read_17/ReadVariableOpRead_17/ReadVariableOp26
Read_18/DisableCopyOnReadRead_18/DisableCopyOnRead20
Read_18/ReadVariableOpRead_18/ReadVariableOp26
Read_19/DisableCopyOnReadRead_19/DisableCopyOnRead20
Read_19/ReadVariableOpRead_19/ReadVariableOp24
Read_2/DisableCopyOnReadRead_2/DisableCopyOnRead2.
Read_2/ReadVariableOpRead_2/ReadVariableOp26
Read_20/DisableCopyOnReadRead_20/DisableCopyOnRead20
Read_20/ReadVariableOpRead_20/ReadVariableOp26
Read_21/DisableCopyOnReadRead_21/DisableCopyOnRead20
Read_21/ReadVariableOpRead_21/ReadVariableOp26
Read_22/DisableCopyOnReadRead_22/DisableCopyOnRead20
Read_22/ReadVariableOpRead_22/ReadVariableOp26
Read_23/DisableCopyOnReadRead_23/DisableCopyOnRead20
Read_23/ReadVariableOpRead_23/ReadVariableOp26
Read_24/DisableCopyOnReadRead_24/DisableCopyOnRead20
Read_24/ReadVariableOpRead_24/ReadVariableOp26
Read_25/DisableCopyOnReadRead_25/DisableCopyOnRead20
Read_25/ReadVariableOpRead_25/ReadVariableOp26
Read_26/DisableCopyOnReadRead_26/DisableCopyOnRead20
Read_26/ReadVariableOpRead_26/ReadVariableOp26
Read_27/DisableCopyOnReadRead_27/DisableCopyOnRead20
Read_27/ReadVariableOpRead_27/ReadVariableOp26
Read_28/DisableCopyOnReadRead_28/DisableCopyOnRead20
Read_28/ReadVariableOpRead_28/ReadVariableOp26
Read_29/DisableCopyOnReadRead_29/DisableCopyOnRead20
Read_29/ReadVariableOpRead_29/ReadVariableOp24
Read_3/DisableCopyOnReadRead_3/DisableCopyOnRead2.
Read_3/ReadVariableOpRead_3/ReadVariableOp24
Read_4/DisableCopyOnReadRead_4/DisableCopyOnRead2.
Read_4/ReadVariableOpRead_4/ReadVariableOp24
Read_5/DisableCopyOnReadRead_5/DisableCopyOnRead2.
Read_5/ReadVariableOpRead_5/ReadVariableOp24
Read_6/DisableCopyOnReadRead_6/DisableCopyOnRead2.
Read_6/ReadVariableOpRead_6/ReadVariableOp24
Read_7/DisableCopyOnReadRead_7/DisableCopyOnRead2.
Read_7/ReadVariableOpRead_7/ReadVariableOp24
Read_8/DisableCopyOnReadRead_8/DisableCopyOnRead2.
Read_8/ReadVariableOpRead_8/ReadVariableOp24
Read_9/DisableCopyOnReadRead_9/DisableCopyOnRead2.
Read_9/ReadVariableOpRead_9/ReadVariableOp:

_output_shapes
: :C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�
e
,__inference_sequential_3_layer_call_fn_11396

inputs
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_3_layer_call_and_return_conditional_losses_10411w
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:���������
'`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
'22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�	
�
0__inference_finetuning_model_layer_call_fn_11235

inputs!
unknown: 
	unknown_0: #
	unknown_1:  
	unknown_2: 
	unknown_3:	�. 
	unknown_4: 
	unknown_5: 

	unknown_6:

identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *T
fORM
K__inference_finetuning_model_layer_call_and_return_conditional_losses_11049o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������
`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:���������
': : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�~
�
!__inference__traced_restore_12160
file_prefix1
assignvariableop_dense_3_kernel: 
-
assignvariableop_1_dense_3_bias:
:
 assignvariableop_2_conv_0_kernel: ,
assignvariableop_3_conv_0_bias: :
 assignvariableop_4_conv_1_kernel:  ,
assignvariableop_5_conv_1_bias: 1
assignvariableop_6_fc_0_kernel:	�. *
assignvariableop_7_fc_0_bias: &
assignvariableop_8_iteration:	 *
 assignvariableop_9_learning_rate: B
(assignvariableop_10_adam_m_conv_0_kernel: B
(assignvariableop_11_adam_v_conv_0_kernel: 4
&assignvariableop_12_adam_m_conv_0_bias: 4
&assignvariableop_13_adam_v_conv_0_bias: B
(assignvariableop_14_adam_m_conv_1_kernel:  B
(assignvariableop_15_adam_v_conv_1_kernel:  4
&assignvariableop_16_adam_m_conv_1_bias: 4
&assignvariableop_17_adam_v_conv_1_bias: 9
&assignvariableop_18_adam_m_fc_0_kernel:	�. 9
&assignvariableop_19_adam_v_fc_0_kernel:	�. 2
$assignvariableop_20_adam_m_fc_0_bias: 2
$assignvariableop_21_adam_v_fc_0_bias: ;
)assignvariableop_22_adam_m_dense_3_kernel: 
;
)assignvariableop_23_adam_v_dense_3_kernel: 
5
'assignvariableop_24_adam_m_dense_3_bias:
5
'assignvariableop_25_adam_v_dense_3_bias:
%
assignvariableop_26_total_1: %
assignvariableop_27_count_1: #
assignvariableop_28_total: #
assignvariableop_29_count: 
identity_31��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_10�AssignVariableOp_11�AssignVariableOp_12�AssignVariableOp_13�AssignVariableOp_14�AssignVariableOp_15�AssignVariableOp_16�AssignVariableOp_17�AssignVariableOp_18�AssignVariableOp_19�AssignVariableOp_2�AssignVariableOp_20�AssignVariableOp_21�AssignVariableOp_22�AssignVariableOp_23�AssignVariableOp_24�AssignVariableOp_25�AssignVariableOp_26�AssignVariableOp_27�AssignVariableOp_28�AssignVariableOp_29�AssignVariableOp_3�AssignVariableOp_4�AssignVariableOp_5�AssignVariableOp_6�AssignVariableOp_7�AssignVariableOp_8�AssignVariableOp_9�
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*�
value�B�B6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB&variables/0/.ATTRIBUTES/VARIABLE_VALUEB&variables/1/.ATTRIBUTES/VARIABLE_VALUEB&variables/2/.ATTRIBUTES/VARIABLE_VALUEB&variables/3/.ATTRIBUTES/VARIABLE_VALUEB&variables/4/.ATTRIBUTES/VARIABLE_VALUEB&variables/5/.ATTRIBUTES/VARIABLE_VALUEB0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUEB3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*Q
valueHBFB B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B �
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*�
_output_shapes~
|:::::::::::::::::::::::::::::::*-
dtypes#
!2	[
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOpAssignVariableOpassignvariableop_dense_3_kernelIdentity:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_1AssignVariableOpassignvariableop_1_dense_3_biasIdentity_1:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_2AssignVariableOp assignvariableop_2_conv_0_kernelIdentity_2:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_3AssignVariableOpassignvariableop_3_conv_0_biasIdentity_3:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_4AssignVariableOp assignvariableop_4_conv_1_kernelIdentity_4:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_5AssignVariableOpassignvariableop_5_conv_1_biasIdentity_5:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_6AssignVariableOpassignvariableop_6_fc_0_kernelIdentity_6:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_7AssignVariableOpassignvariableop_7_fc_0_biasIdentity_7:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0	*
_output_shapes
:�
AssignVariableOp_8AssignVariableOpassignvariableop_8_iterationIdentity_8:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0	]

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_9AssignVariableOp assignvariableop_9_learning_rateIdentity_9:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_10AssignVariableOp(assignvariableop_10_adam_m_conv_0_kernelIdentity_10:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_11AssignVariableOp(assignvariableop_11_adam_v_conv_0_kernelIdentity_11:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_12AssignVariableOp&assignvariableop_12_adam_m_conv_0_biasIdentity_12:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_13AssignVariableOp&assignvariableop_13_adam_v_conv_0_biasIdentity_13:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_14AssignVariableOp(assignvariableop_14_adam_m_conv_1_kernelIdentity_14:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_15AssignVariableOp(assignvariableop_15_adam_v_conv_1_kernelIdentity_15:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_16AssignVariableOp&assignvariableop_16_adam_m_conv_1_biasIdentity_16:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_17AssignVariableOp&assignvariableop_17_adam_v_conv_1_biasIdentity_17:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_18AssignVariableOp&assignvariableop_18_adam_m_fc_0_kernelIdentity_18:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_19AssignVariableOp&assignvariableop_19_adam_v_fc_0_kernelIdentity_19:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_20AssignVariableOp$assignvariableop_20_adam_m_fc_0_biasIdentity_20:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_21AssignVariableOp$assignvariableop_21_adam_v_fc_0_biasIdentity_21:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_22AssignVariableOp)assignvariableop_22_adam_m_dense_3_kernelIdentity_22:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_23AssignVariableOp)assignvariableop_23_adam_v_dense_3_kernelIdentity_23:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_24AssignVariableOp'assignvariableop_24_adam_m_dense_3_biasIdentity_24:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_25AssignVariableOp'assignvariableop_25_adam_v_dense_3_biasIdentity_25:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_26AssignVariableOpassignvariableop_26_total_1Identity_26:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_27AssignVariableOpassignvariableop_27_count_1Identity_27:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_28AssignVariableOpassignvariableop_28_totalIdentity_28:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_29AssignVariableOpassignvariableop_29_countIdentity_29:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0Y
NoOpNoOp"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 �
Identity_30Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: W
Identity_31IdentityIdentity_30:output:0^NoOp_1*
T0*
_output_shapes
: �
NoOp_1NoOp^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9*"
_acd_function_control_output(*
_output_shapes
 "#
identity_31Identity_31:output:0*Q
_input_shapes@
>: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 2*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142*
AssignVariableOp_15AssignVariableOp_152*
AssignVariableOp_16AssignVariableOp_162*
AssignVariableOp_17AssignVariableOp_172*
AssignVariableOp_18AssignVariableOp_182*
AssignVariableOp_19AssignVariableOp_192(
AssignVariableOp_1AssignVariableOp_12*
AssignVariableOp_20AssignVariableOp_202*
AssignVariableOp_21AssignVariableOp_212*
AssignVariableOp_22AssignVariableOp_222*
AssignVariableOp_23AssignVariableOp_232*
AssignVariableOp_24AssignVariableOp_242*
AssignVariableOp_25AssignVariableOp_252*
AssignVariableOp_26AssignVariableOp_262*
AssignVariableOp_27AssignVariableOp_272*
AssignVariableOp_28AssignVariableOp_282*
AssignVariableOp_29AssignVariableOp_292(
AssignVariableOp_2AssignVariableOp_22(
AssignVariableOp_3AssignVariableOp_32(
AssignVariableOp_4AssignVariableOp_42(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_92$
AssignVariableOpAssignVariableOp:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�
�
$__inference_fc_0_layer_call_fn_11784

inputs
unknown:	�. 
	unknown_0: 
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_fc_0_layer_call_and_return_conditional_losses_10539o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:��������� `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������.: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:����������.
 
_user_specified_nameinputs
�
^
B__inference_flatten_layer_call_and_return_conditional_losses_11775

inputs
identityV
ConstConst*
_output_shapes
:*
dtype0*
valueB"����@  ]
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:����������.Y
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:����������."
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:��������� :W S
/
_output_shapes
:��������� 
 
_user_specified_nameinputs
�
e
I__inference_gaussian_noise_layer_call_and_return_conditional_losses_10576

inputs
identityV
IdentityIdentityinputs*
T0*/
_output_shapes
:���������
'"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
':W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�	
�
0__inference_finetuning_model_layer_call_fn_11000
input_5!
unknown: 
	unknown_0: #
	unknown_1:  
	unknown_2: 
	unknown_3:	�. 
	unknown_4: 
	unknown_5: 

	unknown_6:

identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_5unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *T
fORM
K__inference_finetuning_model_layer_call_and_return_conditional_losses_10981o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������
`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:���������
': : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:X T
/
_output_shapes
:���������
'
!
_user_specified_name	input_5
�
�
A__inference_conv_0_layer_call_and_return_conditional_losses_10481

inputs8
conv2d_readvariableop_resource: -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�Conv2D/ReadVariableOp�-conv_0/bias/Regularizer/L2Loss/ReadVariableOp�/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp|
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
: *
dtype0�
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������# *
paddingVALID*
strides
r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0}
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������# X
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:���������# �
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
: *
dtype0�
 conv_0/kernel/Regularizer/L2LossL2Loss7conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/kernel/Regularizer/mulMul(conv_0/kernel/Regularizer/mul/x:output:0)conv_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-conv_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_0/bias/Regularizer/L2LossL2Loss5conv_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/bias/Regularizer/mulMul&conv_0/bias/Regularizer/mul/x:output:0'conv_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: i
IdentityIdentityRelu:activations:0^NoOp*
T0*/
_output_shapes
:���������# �
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp.^conv_0/bias/Regularizer/L2Loss/ReadVariableOp0^conv_0/kernel/Regularizer/L2Loss/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������
': : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp2^
-conv_0/bias/Regularizer/L2Loss/ReadVariableOp-conv_0/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�
^
B__inference_flatten_layer_call_and_return_conditional_losses_10518

inputs
identityV
ConstConst*
_output_shapes
:*
dtype0*
valueB"����@  ]
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:����������.Y
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:����������."
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:��������� :W S
/
_output_shapes
:��������� 
 
_user_specified_nameinputs
�8
�
G__inference_sequential_2_layer_call_and_return_conditional_losses_10729

inputs&
conv_0_10688: 
conv_0_10690: &
conv_1_10693:  
conv_1_10695: 

fc_0_10699:	�. 

fc_0_10701: 
identity��conv_0/StatefulPartitionedCall�-conv_0/bias/Regularizer/L2Loss/ReadVariableOp�/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp�conv_1/StatefulPartitionedCall�-conv_1/bias/Regularizer/L2Loss/ReadVariableOp�/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp�fc_0/StatefulPartitionedCall�+fc_0/bias/Regularizer/L2Loss/ReadVariableOp�-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp�
gaussian_noise/PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_gaussian_noise_layer_call_and_return_conditional_losses_10576�
conv_0/StatefulPartitionedCallStatefulPartitionedCall'gaussian_noise/PartitionedCall:output:0conv_0_10688conv_0_10690*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������# *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *J
fERC
A__inference_conv_0_layer_call_and_return_conditional_losses_10481�
conv_1/StatefulPartitionedCallStatefulPartitionedCall'conv_0/StatefulPartitionedCall:output:0conv_1_10693conv_1_10695*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *J
fERC
A__inference_conv_1_layer_call_and_return_conditional_losses_10506�
flatten/PartitionedCallPartitionedCall'conv_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������.* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_flatten_layer_call_and_return_conditional_losses_10518�
fc_0/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0
fc_0_10699
fc_0_10701*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_fc_0_layer_call_and_return_conditional_losses_10539�
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv_0_10688*&
_output_shapes
: *
dtype0�
 conv_0/kernel/Regularizer/L2LossL2Loss7conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/kernel/Regularizer/mulMul(conv_0/kernel/Regularizer/mul/x:output:0)conv_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: v
-conv_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv_0_10690*
_output_shapes
: *
dtype0�
conv_0/bias/Regularizer/L2LossL2Loss5conv_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/bias/Regularizer/mulMul&conv_0/bias/Regularizer/mul/x:output:0'conv_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv_1_10693*&
_output_shapes
:  *
dtype0�
 conv_1/kernel/Regularizer/L2LossL2Loss7conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/kernel/Regularizer/mulMul(conv_1/kernel/Regularizer/mul/x:output:0)conv_1/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: v
-conv_1/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv_1_10695*
_output_shapes
: *
dtype0�
conv_1/bias/Regularizer/L2LossL2Loss5conv_1/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/bias/Regularizer/mulMul&conv_1/bias/Regularizer/mul/x:output:0'conv_1/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: y
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp
fc_0_10699*
_output_shapes
:	�. *
dtype0�
fc_0/kernel/Regularizer/L2LossL2Loss5fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
fc_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/kernel/Regularizer/mulMul&fc_0/kernel/Regularizer/mul/x:output:0'fc_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: r
+fc_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp
fc_0_10701*
_output_shapes
: *
dtype0|
fc_0/bias/Regularizer/L2LossL2Loss3fc_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: `
fc_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/bias/Regularizer/mulMul$fc_0/bias/Regularizer/mul/x:output:0%fc_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: t
IdentityIdentity%fc_0/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:��������� �
NoOpNoOp^conv_0/StatefulPartitionedCall.^conv_0/bias/Regularizer/L2Loss/ReadVariableOp0^conv_0/kernel/Regularizer/L2Loss/ReadVariableOp^conv_1/StatefulPartitionedCall.^conv_1/bias/Regularizer/L2Loss/ReadVariableOp0^conv_1/kernel/Regularizer/L2Loss/ReadVariableOp^fc_0/StatefulPartitionedCall,^fc_0/bias/Regularizer/L2Loss/ReadVariableOp.^fc_0/kernel/Regularizer/L2Loss/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*:
_input_shapes)
':���������
': : : : : : 2@
conv_0/StatefulPartitionedCallconv_0/StatefulPartitionedCall2^
-conv_0/bias/Regularizer/L2Loss/ReadVariableOp-conv_0/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp2@
conv_1/StatefulPartitionedCallconv_1/StatefulPartitionedCall2^
-conv_1/bias/Regularizer/L2Loss/ReadVariableOp-conv_1/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp2<
fc_0/StatefulPartitionedCallfc_0/StatefulPartitionedCall2Z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOp+fc_0/bias/Regularizer/L2Loss/ReadVariableOp2^
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�	
�
B__inference_dense_3_layer_call_and_return_conditional_losses_10853

inputs0
matmul_readvariableop_resource: 
-
biasadd_readvariableop_resource:

identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

: 
*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:
*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������
w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:��������� 
 
_user_specified_nameinputs
�
J
.__inference_gaussian_noise_layer_call_fn_11693

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_gaussian_noise_layer_call_and_return_conditional_losses_10576h
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:���������
'"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
':W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�5
�
K__inference_finetuning_model_layer_call_and_return_conditional_losses_10884
input_5,
sequential_2_10830:  
sequential_2_10832: ,
sequential_2_10834:   
sequential_2_10836: %
sequential_2_10838:	�.  
sequential_2_10840: 
dense_3_10854: 

dense_3_10856:

identity��-conv_0/bias/Regularizer/L2Loss/ReadVariableOp�/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp�-conv_1/bias/Regularizer/L2Loss/ReadVariableOp�/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp�dense_3/StatefulPartitionedCall�+fc_0/bias/Regularizer/L2Loss/ReadVariableOp�-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp�$sequential_2/StatefulPartitionedCall�$sequential_3/StatefulPartitionedCall�
$sequential_3/StatefulPartitionedCallStatefulPartitionedCallinput_5*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_3_layer_call_and_return_conditional_losses_10411�
$sequential_2/StatefulPartitionedCallStatefulPartitionedCall-sequential_3/StatefulPartitionedCall:output:0sequential_2_10830sequential_2_10832sequential_2_10834sequential_2_10836sequential_2_10838sequential_2_10840*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *(
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_2_layer_call_and_return_conditional_losses_10667�
dense_3/StatefulPartitionedCallStatefulPartitionedCall-sequential_2/StatefulPartitionedCall:output:0dense_3_10854dense_3_10856*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_3_layer_call_and_return_conditional_losses_10853�
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10830*&
_output_shapes
: *
dtype0�
 conv_0/kernel/Regularizer/L2LossL2Loss7conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/kernel/Regularizer/mulMul(conv_0/kernel/Regularizer/mul/x:output:0)conv_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: |
-conv_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10832*
_output_shapes
: *
dtype0�
conv_0/bias/Regularizer/L2LossL2Loss5conv_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/bias/Regularizer/mulMul&conv_0/bias/Regularizer/mul/x:output:0'conv_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10834*&
_output_shapes
:  *
dtype0�
 conv_1/kernel/Regularizer/L2LossL2Loss7conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/kernel/Regularizer/mulMul(conv_1/kernel/Regularizer/mul/x:output:0)conv_1/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: |
-conv_1/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10836*
_output_shapes
: *
dtype0�
conv_1/bias/Regularizer/L2LossL2Loss5conv_1/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/bias/Regularizer/mulMul&conv_1/bias/Regularizer/mul/x:output:0'conv_1/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10838*
_output_shapes
:	�. *
dtype0�
fc_0/kernel/Regularizer/L2LossL2Loss5fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
fc_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/kernel/Regularizer/mulMul&fc_0/kernel/Regularizer/mul/x:output:0'fc_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10840*
_output_shapes
: *
dtype0|
fc_0/bias/Regularizer/L2LossL2Loss3fc_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: `
fc_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/bias/Regularizer/mulMul$fc_0/bias/Regularizer/mul/x:output:0%fc_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: w
IdentityIdentity(dense_3/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������
�
NoOpNoOp.^conv_0/bias/Regularizer/L2Loss/ReadVariableOp0^conv_0/kernel/Regularizer/L2Loss/ReadVariableOp.^conv_1/bias/Regularizer/L2Loss/ReadVariableOp0^conv_1/kernel/Regularizer/L2Loss/ReadVariableOp ^dense_3/StatefulPartitionedCall,^fc_0/bias/Regularizer/L2Loss/ReadVariableOp.^fc_0/kernel/Regularizer/L2Loss/ReadVariableOp%^sequential_2/StatefulPartitionedCall%^sequential_3/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:���������
': : : : : : : : 2^
-conv_0/bias/Regularizer/L2Loss/ReadVariableOp-conv_0/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp2^
-conv_1/bias/Regularizer/L2Loss/ReadVariableOp-conv_1/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2Z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOp+fc_0/bias/Regularizer/L2Loss/ReadVariableOp2^
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp2L
$sequential_2/StatefulPartitionedCall$sequential_2/StatefulPartitionedCall2L
$sequential_3/StatefulPartitionedCall$sequential_3/StatefulPartitionedCall:X T
/
_output_shapes
:���������
'
!
_user_specified_name	input_5
�
n
5__inference_random_color_affine_2_layer_call_fn_11635

images
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallimages*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *Y
fTRR
P__inference_random_color_affine_2_layer_call_and_return_conditional_losses_10391w
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:���������
'`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
'22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������
'
 
_user_specified_nameimages
�	
�
B__inference_dense_3_layer_call_and_return_conditional_losses_11630

inputs0
matmul_readvariableop_resource: 
-
biasadd_readvariableop_resource:

identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

: 
*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:
*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
_
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������
w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:��������� 
 
_user_specified_nameinputs
�
C
'__inference_flatten_layer_call_fn_11769

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������.* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_flatten_layer_call_and_return_conditional_losses_10518a
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:����������."
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:��������� :W S
/
_output_shapes
:��������� 
 
_user_specified_nameinputs
�
e
I__inference_gaussian_noise_layer_call_and_return_conditional_losses_11708

inputs
identityV
IdentityIdentityinputs*
T0*/
_output_shapes
:���������
'"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
':W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�@
�
G__inference_sequential_2_layer_call_and_return_conditional_losses_11611

inputs?
%conv_0_conv2d_readvariableop_resource: 4
&conv_0_biasadd_readvariableop_resource: ?
%conv_1_conv2d_readvariableop_resource:  4
&conv_1_biasadd_readvariableop_resource: 6
#fc_0_matmul_readvariableop_resource:	�. 2
$fc_0_biasadd_readvariableop_resource: 
identity��conv_0/BiasAdd/ReadVariableOp�conv_0/Conv2D/ReadVariableOp�-conv_0/bias/Regularizer/L2Loss/ReadVariableOp�/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp�conv_1/BiasAdd/ReadVariableOp�conv_1/Conv2D/ReadVariableOp�-conv_1/bias/Regularizer/L2Loss/ReadVariableOp�/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp�fc_0/BiasAdd/ReadVariableOp�fc_0/MatMul/ReadVariableOp�+fc_0/bias/Regularizer/L2Loss/ReadVariableOp�-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp�
conv_0/Conv2D/ReadVariableOpReadVariableOp%conv_0_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype0�
conv_0/Conv2DConv2Dinputs$conv_0/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������# *
paddingVALID*
strides
�
conv_0/BiasAdd/ReadVariableOpReadVariableOp&conv_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_0/BiasAddBiasAddconv_0/Conv2D:output:0%conv_0/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������# f
conv_0/ReluReluconv_0/BiasAdd:output:0*
T0*/
_output_shapes
:���������# �
conv_1/Conv2D/ReadVariableOpReadVariableOp%conv_1_conv2d_readvariableop_resource*&
_output_shapes
:  *
dtype0�
conv_1/Conv2DConv2Dconv_0/Relu:activations:0$conv_1/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:��������� *
paddingVALID*
strides
�
conv_1/BiasAdd/ReadVariableOpReadVariableOp&conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_1/BiasAddBiasAddconv_1/Conv2D:output:0%conv_1/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:��������� f
conv_1/ReluReluconv_1/BiasAdd:output:0*
T0*/
_output_shapes
:��������� ^
flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"����@  �
flatten/ReshapeReshapeconv_1/Relu:activations:0flatten/Const:output:0*
T0*(
_output_shapes
:����������.
fc_0/MatMul/ReadVariableOpReadVariableOp#fc_0_matmul_readvariableop_resource*
_output_shapes
:	�. *
dtype0�
fc_0/MatMulMatMulflatten/Reshape:output:0"fc_0/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� |
fc_0/BiasAdd/ReadVariableOpReadVariableOp$fc_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
fc_0/BiasAddBiasAddfc_0/MatMul:product:0#fc_0/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� Z
	fc_0/ReluRelufc_0/BiasAdd:output:0*
T0*'
_output_shapes
:��������� �
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp%conv_0_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype0�
 conv_0/kernel/Regularizer/L2LossL2Loss7conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/kernel/Regularizer/mulMul(conv_0/kernel/Regularizer/mul/x:output:0)conv_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-conv_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp&conv_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_0/bias/Regularizer/L2LossL2Loss5conv_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/bias/Regularizer/mulMul&conv_0/bias/Regularizer/mul/x:output:0'conv_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp%conv_1_conv2d_readvariableop_resource*&
_output_shapes
:  *
dtype0�
 conv_1/kernel/Regularizer/L2LossL2Loss7conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/kernel/Regularizer/mulMul(conv_1/kernel/Regularizer/mul/x:output:0)conv_1/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-conv_1/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp&conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_1/bias/Regularizer/L2LossL2Loss5conv_1/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/bias/Regularizer/mulMul&conv_1/bias/Regularizer/mul/x:output:0'conv_1/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp#fc_0_matmul_readvariableop_resource*
_output_shapes
:	�. *
dtype0�
fc_0/kernel/Regularizer/L2LossL2Loss5fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
fc_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/kernel/Regularizer/mulMul&fc_0/kernel/Regularizer/mul/x:output:0'fc_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
+fc_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp$fc_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0|
fc_0/bias/Regularizer/L2LossL2Loss3fc_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: `
fc_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/bias/Regularizer/mulMul$fc_0/bias/Regularizer/mul/x:output:0%fc_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: f
IdentityIdentityfc_0/Relu:activations:0^NoOp*
T0*'
_output_shapes
:��������� �
NoOpNoOp^conv_0/BiasAdd/ReadVariableOp^conv_0/Conv2D/ReadVariableOp.^conv_0/bias/Regularizer/L2Loss/ReadVariableOp0^conv_0/kernel/Regularizer/L2Loss/ReadVariableOp^conv_1/BiasAdd/ReadVariableOp^conv_1/Conv2D/ReadVariableOp.^conv_1/bias/Regularizer/L2Loss/ReadVariableOp0^conv_1/kernel/Regularizer/L2Loss/ReadVariableOp^fc_0/BiasAdd/ReadVariableOp^fc_0/MatMul/ReadVariableOp,^fc_0/bias/Regularizer/L2Loss/ReadVariableOp.^fc_0/kernel/Regularizer/L2Loss/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*:
_input_shapes)
':���������
': : : : : : 2>
conv_0/BiasAdd/ReadVariableOpconv_0/BiasAdd/ReadVariableOp2<
conv_0/Conv2D/ReadVariableOpconv_0/Conv2D/ReadVariableOp2^
-conv_0/bias/Regularizer/L2Loss/ReadVariableOp-conv_0/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp2>
conv_1/BiasAdd/ReadVariableOpconv_1/BiasAdd/ReadVariableOp2<
conv_1/Conv2D/ReadVariableOpconv_1/Conv2D/ReadVariableOp2^
-conv_1/bias/Regularizer/L2Loss/ReadVariableOp-conv_1/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp2:
fc_0/BiasAdd/ReadVariableOpfc_0/BiasAdd/ReadVariableOp28
fc_0/MatMul/ReadVariableOpfc_0/MatMul/ReadVariableOp2Z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOp+fc_0/bias/Regularizer/L2Loss/ReadVariableOp2^
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�
f
,__inference_sequential_3_layer_call_fn_10414
input_6
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_6*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_3_layer_call_and_return_conditional_losses_10411w
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:���������
'`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
'22
StatefulPartitionedCallStatefulPartitionedCall:X T
/
_output_shapes
:���������
'
!
_user_specified_name	input_6
�
I
,__inference_sequential_3_layer_call_fn_10424
input_6
identity�
PartitionedCallPartitionedCallinput_6*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_3_layer_call_and_return_conditional_losses_10421h
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:���������
'"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
':X T
/
_output_shapes
:���������
'
!
_user_specified_name	input_6
�3
�
K__inference_finetuning_model_layer_call_and_return_conditional_losses_10931
input_5,
sequential_2_10888:  
sequential_2_10890: ,
sequential_2_10892:   
sequential_2_10894: %
sequential_2_10896:	�.  
sequential_2_10898: 
dense_3_10901: 

dense_3_10903:

identity��-conv_0/bias/Regularizer/L2Loss/ReadVariableOp�/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp�-conv_1/bias/Regularizer/L2Loss/ReadVariableOp�/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp�dense_3/StatefulPartitionedCall�+fc_0/bias/Regularizer/L2Loss/ReadVariableOp�-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp�$sequential_2/StatefulPartitionedCall�
sequential_3/PartitionedCallPartitionedCallinput_5*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_3_layer_call_and_return_conditional_losses_10421�
$sequential_2/StatefulPartitionedCallStatefulPartitionedCall%sequential_3/PartitionedCall:output:0sequential_2_10888sequential_2_10890sequential_2_10892sequential_2_10894sequential_2_10896sequential_2_10898*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *(
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_2_layer_call_and_return_conditional_losses_10729�
dense_3/StatefulPartitionedCallStatefulPartitionedCall-sequential_2/StatefulPartitionedCall:output:0dense_3_10901dense_3_10903*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_3_layer_call_and_return_conditional_losses_10853�
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10888*&
_output_shapes
: *
dtype0�
 conv_0/kernel/Regularizer/L2LossL2Loss7conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/kernel/Regularizer/mulMul(conv_0/kernel/Regularizer/mul/x:output:0)conv_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: |
-conv_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10890*
_output_shapes
: *
dtype0�
conv_0/bias/Regularizer/L2LossL2Loss5conv_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/bias/Regularizer/mulMul&conv_0/bias/Regularizer/mul/x:output:0'conv_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10892*&
_output_shapes
:  *
dtype0�
 conv_1/kernel/Regularizer/L2LossL2Loss7conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/kernel/Regularizer/mulMul(conv_1/kernel/Regularizer/mul/x:output:0)conv_1/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: |
-conv_1/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10894*
_output_shapes
: *
dtype0�
conv_1/bias/Regularizer/L2LossL2Loss5conv_1/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/bias/Regularizer/mulMul&conv_1/bias/Regularizer/mul/x:output:0'conv_1/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10896*
_output_shapes
:	�. *
dtype0�
fc_0/kernel/Regularizer/L2LossL2Loss5fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
fc_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/kernel/Regularizer/mulMul&fc_0/kernel/Regularizer/mul/x:output:0'fc_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_10898*
_output_shapes
: *
dtype0|
fc_0/bias/Regularizer/L2LossL2Loss3fc_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: `
fc_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/bias/Regularizer/mulMul$fc_0/bias/Regularizer/mul/x:output:0%fc_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: w
IdentityIdentity(dense_3/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������
�
NoOpNoOp.^conv_0/bias/Regularizer/L2Loss/ReadVariableOp0^conv_0/kernel/Regularizer/L2Loss/ReadVariableOp.^conv_1/bias/Regularizer/L2Loss/ReadVariableOp0^conv_1/kernel/Regularizer/L2Loss/ReadVariableOp ^dense_3/StatefulPartitionedCall,^fc_0/bias/Regularizer/L2Loss/ReadVariableOp.^fc_0/kernel/Regularizer/L2Loss/ReadVariableOp%^sequential_2/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:���������
': : : : : : : : 2^
-conv_0/bias/Regularizer/L2Loss/ReadVariableOp-conv_0/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp2^
-conv_1/bias/Regularizer/L2Loss/ReadVariableOp-conv_1/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2Z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOp+fc_0/bias/Regularizer/L2Loss/ReadVariableOp2^
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp2L
$sequential_2/StatefulPartitionedCall$sequential_2/StatefulPartitionedCall:X T
/
_output_shapes
:���������
'
!
_user_specified_name	input_5
�	
h
I__inference_gaussian_noise_layer_call_and_return_conditional_losses_10460

inputs
identity�I
ShapeShapeinputs*
T0*
_output_shapes
::��W
random_normal/meanConst*
_output_shapes
: *
dtype0*
valueB
 *    Y
random_normal/stddevConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
"random_normal/RandomStandardNormalRandomStandardNormalShape:output:0*
T0*/
_output_shapes
:���������
'*
dtype0�
random_normal/mulMul+random_normal/RandomStandardNormal:output:0random_normal/stddev:output:0*
T0*/
_output_shapes
:���������
'�
random_normalAddV2random_normal/mul:z:0random_normal/mean:output:0*
T0*/
_output_shapes
:���������
'a
addAddV2inputsrandom_normal:z:0*
T0*/
_output_shapes
:���������
'W
IdentityIdentityadd:z:0*
T0*/
_output_shapes
:���������
'"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
':W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�#
o
P__inference_random_color_affine_2_layer_call_and_return_conditional_losses_11679

images
identity�I
ShapeShapeimages*
T0*
_output_shapes
::��]
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: _
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:_
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:�
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_maskX
random_uniform/shape/1Const*
_output_shapes
: *
dtype0*
value	B :X
random_uniform/shape/2Const*
_output_shapes
: *
dtype0*
value	B :X
random_uniform/shape/3Const*
_output_shapes
: *
dtype0*
value	B :�
random_uniform/shapePackstrided_slice:output:0random_uniform/shape/1:output:0random_uniform/shape/2:output:0random_uniform/shape/3:output:0*
N*
T0*
_output_shapes
:W
random_uniform/minConst*
_output_shapes
: *
dtype0*
valueB
 *����W
random_uniform/maxConst*
_output_shapes
: *
dtype0*
valueB
 *���>�
random_uniform/RandomUniformRandomUniformrandom_uniform/shape:output:0*
T0*/
_output_shapes
:���������*
dtype0t
random_uniform/subSubrandom_uniform/max:output:0random_uniform/min:output:0*
T0*
_output_shapes
: �
random_uniform/mulMul%random_uniform/RandomUniform:output:0random_uniform/sub:z:0*
T0*/
_output_shapes
:����������
random_uniformAddV2random_uniform/mul:z:0random_uniform/min:output:0*
T0*/
_output_shapes
:���������b
AddAddV2imagesrandom_uniform:z:0*
T0*/
_output_shapes
:���������
'Z
random_uniform_1/shape/1Const*
_output_shapes
: *
dtype0*
value	B :Z
random_uniform_1/shape/2Const*
_output_shapes
: *
dtype0*
value	B :Z
random_uniform_1/shape/3Const*
_output_shapes
: *
dtype0*
value	B :�
random_uniform_1/shapePackstrided_slice:output:0!random_uniform_1/shape/1:output:0!random_uniform_1/shape/2:output:0!random_uniform_1/shape/3:output:0*
N*
T0*
_output_shapes
:Y
random_uniform_1/minConst*
_output_shapes
: *
dtype0*
valueB
 *��̽Y
random_uniform_1/maxConst*
_output_shapes
: *
dtype0*
valueB
 *���=�
random_uniform_1/RandomUniformRandomUniformrandom_uniform_1/shape:output:0*
T0*/
_output_shapes
:���������*
dtype0z
random_uniform_1/subSubrandom_uniform_1/max:output:0random_uniform_1/min:output:0*
T0*
_output_shapes
: �
random_uniform_1/mulMul'random_uniform_1/RandomUniform:output:0random_uniform_1/sub:z:0*
T0*/
_output_shapes
:����������
random_uniform_1AddV2random_uniform_1/mul:z:0random_uniform_1/min:output:0*
T0*/
_output_shapes
:���������k
Mean/reduction_indicesConst*
_output_shapes
:*
dtype0*!
valueB"         �
MeanMeanAdd:z:0Mean/reduction_indices:output:0*
T0*/
_output_shapes
:���������*
	keep_dims(\
subSubAdd:z:0Mean:output:0*
T0*/
_output_shapes
:���������
'c
mulMulrandom_uniform_1:z:0sub:z:0*
T0*/
_output_shapes
:���������
'Z
Add_1AddV2Add:z:0mul:z:0*
T0*/
_output_shapes
:���������
'\
clip_by_value/Minimum/yConst*
_output_shapes
: *
dtype0*
valueB
 *  �?�
clip_by_value/MinimumMinimum	Add_1:z:0 clip_by_value/Minimum/y:output:0*
T0*/
_output_shapes
:���������
'T
clip_by_value/yConst*
_output_shapes
: *
dtype0*
valueB
 *    �
clip_by_valueMaximumclip_by_value/Minimum:z:0clip_by_value/y:output:0*
T0*/
_output_shapes
:���������
'a
IdentityIdentityclip_by_value:z:0*
T0*/
_output_shapes
:���������
'"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
':W S
/
_output_shapes
:���������
'
 
_user_specified_nameimages
�
�
A__inference_conv_1_layer_call_and_return_conditional_losses_10506

inputs8
conv2d_readvariableop_resource:  -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�Conv2D/ReadVariableOp�-conv_1/bias/Regularizer/L2Loss/ReadVariableOp�/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp|
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:  *
dtype0�
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:��������� *
paddingVALID*
strides
r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0}
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:��������� X
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:��������� �
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:  *
dtype0�
 conv_1/kernel/Regularizer/L2LossL2Loss7conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/kernel/Regularizer/mulMul(conv_1/kernel/Regularizer/mul/x:output:0)conv_1/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-conv_1/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_1/bias/Regularizer/L2LossL2Loss5conv_1/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/bias/Regularizer/mulMul&conv_1/bias/Regularizer/mul/x:output:0'conv_1/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: i
IdentityIdentityRelu:activations:0^NoOp*
T0*/
_output_shapes
:��������� �
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp.^conv_1/bias/Regularizer/L2Loss/ReadVariableOp0^conv_1/kernel/Regularizer/L2Loss/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������# : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp2^
-conv_1/bias/Regularizer/L2Loss/ReadVariableOp-conv_1/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:W S
/
_output_shapes
:���������# 
 
_user_specified_nameinputs
�
c
G__inference_sequential_3_layer_call_and_return_conditional_losses_10421

inputs
identity�
%random_color_affine_2/PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *Y
fTRR
P__inference_random_color_affine_2_layer_call_and_return_conditional_losses_10400~
IdentityIdentity.random_color_affine_2/PartitionedCall:output:0*
T0*/
_output_shapes
:���������
'"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
':W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�
�
__inference_loss_fn_4_11848I
6fc_0_kernel_regularizer_l2loss_readvariableop_resource:	�. 
identity��-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp�
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp6fc_0_kernel_regularizer_l2loss_readvariableop_resource*
_output_shapes
:	�. *
dtype0�
fc_0/kernel/Regularizer/L2LossL2Loss5fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
fc_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/kernel/Regularizer/mulMul&fc_0/kernel/Regularizer/mul/x:output:0'fc_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: ]
IdentityIdentityfc_0/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: v
NoOpNoOp.^fc_0/kernel/Regularizer/L2Loss/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2^
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp
�
�
&__inference_conv_0_layer_call_fn_11717

inputs!
unknown: 
	unknown_0: 
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������# *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *J
fERC
A__inference_conv_0_layer_call_and_return_conditional_losses_10481w
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:���������# `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������
': : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�9
�
G__inference_sequential_2_layer_call_and_return_conditional_losses_10570
gaussian_noise_input&
conv_0_10482: 
conv_0_10484: &
conv_1_10507:  
conv_1_10509: 

fc_0_10540:	�. 

fc_0_10542: 
identity��conv_0/StatefulPartitionedCall�-conv_0/bias/Regularizer/L2Loss/ReadVariableOp�/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp�conv_1/StatefulPartitionedCall�-conv_1/bias/Regularizer/L2Loss/ReadVariableOp�/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp�fc_0/StatefulPartitionedCall�+fc_0/bias/Regularizer/L2Loss/ReadVariableOp�-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp�&gaussian_noise/StatefulPartitionedCall�
&gaussian_noise/StatefulPartitionedCallStatefulPartitionedCallgaussian_noise_input*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_gaussian_noise_layer_call_and_return_conditional_losses_10460�
conv_0/StatefulPartitionedCallStatefulPartitionedCall/gaussian_noise/StatefulPartitionedCall:output:0conv_0_10482conv_0_10484*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������# *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *J
fERC
A__inference_conv_0_layer_call_and_return_conditional_losses_10481�
conv_1/StatefulPartitionedCallStatefulPartitionedCall'conv_0/StatefulPartitionedCall:output:0conv_1_10507conv_1_10509*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *J
fERC
A__inference_conv_1_layer_call_and_return_conditional_losses_10506�
flatten/PartitionedCallPartitionedCall'conv_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������.* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_flatten_layer_call_and_return_conditional_losses_10518�
fc_0/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0
fc_0_10540
fc_0_10542*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_fc_0_layer_call_and_return_conditional_losses_10539�
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv_0_10482*&
_output_shapes
: *
dtype0�
 conv_0/kernel/Regularizer/L2LossL2Loss7conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/kernel/Regularizer/mulMul(conv_0/kernel/Regularizer/mul/x:output:0)conv_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: v
-conv_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv_0_10484*
_output_shapes
: *
dtype0�
conv_0/bias/Regularizer/L2LossL2Loss5conv_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/bias/Regularizer/mulMul&conv_0/bias/Regularizer/mul/x:output:0'conv_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv_1_10507*&
_output_shapes
:  *
dtype0�
 conv_1/kernel/Regularizer/L2LossL2Loss7conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/kernel/Regularizer/mulMul(conv_1/kernel/Regularizer/mul/x:output:0)conv_1/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: v
-conv_1/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv_1_10509*
_output_shapes
: *
dtype0�
conv_1/bias/Regularizer/L2LossL2Loss5conv_1/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/bias/Regularizer/mulMul&conv_1/bias/Regularizer/mul/x:output:0'conv_1/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: y
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp
fc_0_10540*
_output_shapes
:	�. *
dtype0�
fc_0/kernel/Regularizer/L2LossL2Loss5fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
fc_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/kernel/Regularizer/mulMul&fc_0/kernel/Regularizer/mul/x:output:0'fc_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: r
+fc_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp
fc_0_10542*
_output_shapes
: *
dtype0|
fc_0/bias/Regularizer/L2LossL2Loss3fc_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: `
fc_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/bias/Regularizer/mulMul$fc_0/bias/Regularizer/mul/x:output:0%fc_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: t
IdentityIdentity%fc_0/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:��������� �
NoOpNoOp^conv_0/StatefulPartitionedCall.^conv_0/bias/Regularizer/L2Loss/ReadVariableOp0^conv_0/kernel/Regularizer/L2Loss/ReadVariableOp^conv_1/StatefulPartitionedCall.^conv_1/bias/Regularizer/L2Loss/ReadVariableOp0^conv_1/kernel/Regularizer/L2Loss/ReadVariableOp^fc_0/StatefulPartitionedCall,^fc_0/bias/Regularizer/L2Loss/ReadVariableOp.^fc_0/kernel/Regularizer/L2Loss/ReadVariableOp'^gaussian_noise/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*:
_input_shapes)
':���������
': : : : : : 2@
conv_0/StatefulPartitionedCallconv_0/StatefulPartitionedCall2^
-conv_0/bias/Regularizer/L2Loss/ReadVariableOp-conv_0/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp2@
conv_1/StatefulPartitionedCallconv_1/StatefulPartitionedCall2^
-conv_1/bias/Regularizer/L2Loss/ReadVariableOp-conv_1/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp2<
fc_0/StatefulPartitionedCallfc_0/StatefulPartitionedCall2Z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOp+fc_0/bias/Regularizer/L2Loss/ReadVariableOp2^
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp2P
&gaussian_noise/StatefulPartitionedCall&gaussian_noise/StatefulPartitionedCall:e a
/
_output_shapes
:���������
'
.
_user_specified_namegaussian_noise_input
�
�
?__inference_fc_0_layer_call_and_return_conditional_losses_10539

inputs1
matmul_readvariableop_resource:	�. -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�+fc_0/bias/Regularizer/L2Loss/ReadVariableOp�-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�. *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:��������� �
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�. *
dtype0�
fc_0/kernel/Regularizer/L2LossL2Loss5fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
fc_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/kernel/Regularizer/mulMul&fc_0/kernel/Regularizer/mul/x:output:0'fc_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
+fc_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0|
fc_0/bias/Regularizer/L2LossL2Loss3fc_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: `
fc_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/bias/Regularizer/mulMul$fc_0/bias/Regularizer/mul/x:output:0%fc_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:��������� �
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp,^fc_0/bias/Regularizer/L2Loss/ReadVariableOp.^fc_0/kernel/Regularizer/L2Loss/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������.: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2Z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOp+fc_0/bias/Regularizer/L2Loss/ReadVariableOp2^
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:P L
(
_output_shapes
:����������.
 
_user_specified_nameinputs
�
�
G__inference_sequential_3_layer_call_and_return_conditional_losses_10411

inputs
identity��-random_color_affine_2/StatefulPartitionedCall�
-random_color_affine_2/StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *Y
fTRR
P__inference_random_color_affine_2_layer_call_and_return_conditional_losses_10391�
IdentityIdentity6random_color_affine_2/StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:���������
'v
NoOpNoOp.^random_color_affine_2/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
'2^
-random_color_affine_2/StatefulPartitionedCall-random_color_affine_2/StatefulPartitionedCall:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�9
�
G__inference_sequential_2_layer_call_and_return_conditional_losses_10667

inputs&
conv_0_10626: 
conv_0_10628: &
conv_1_10631:  
conv_1_10633: 

fc_0_10637:	�. 

fc_0_10639: 
identity��conv_0/StatefulPartitionedCall�-conv_0/bias/Regularizer/L2Loss/ReadVariableOp�/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp�conv_1/StatefulPartitionedCall�-conv_1/bias/Regularizer/L2Loss/ReadVariableOp�/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp�fc_0/StatefulPartitionedCall�+fc_0/bias/Regularizer/L2Loss/ReadVariableOp�-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp�&gaussian_noise/StatefulPartitionedCall�
&gaussian_noise/StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_gaussian_noise_layer_call_and_return_conditional_losses_10460�
conv_0/StatefulPartitionedCallStatefulPartitionedCall/gaussian_noise/StatefulPartitionedCall:output:0conv_0_10626conv_0_10628*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������# *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *J
fERC
A__inference_conv_0_layer_call_and_return_conditional_losses_10481�
conv_1/StatefulPartitionedCallStatefulPartitionedCall'conv_0/StatefulPartitionedCall:output:0conv_1_10631conv_1_10633*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *J
fERC
A__inference_conv_1_layer_call_and_return_conditional_losses_10506�
flatten/PartitionedCallPartitionedCall'conv_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������.* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_flatten_layer_call_and_return_conditional_losses_10518�
fc_0/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0
fc_0_10637
fc_0_10639*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_fc_0_layer_call_and_return_conditional_losses_10539�
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv_0_10626*&
_output_shapes
: *
dtype0�
 conv_0/kernel/Regularizer/L2LossL2Loss7conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/kernel/Regularizer/mulMul(conv_0/kernel/Regularizer/mul/x:output:0)conv_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: v
-conv_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv_0_10628*
_output_shapes
: *
dtype0�
conv_0/bias/Regularizer/L2LossL2Loss5conv_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/bias/Regularizer/mulMul&conv_0/bias/Regularizer/mul/x:output:0'conv_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv_1_10631*&
_output_shapes
:  *
dtype0�
 conv_1/kernel/Regularizer/L2LossL2Loss7conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/kernel/Regularizer/mulMul(conv_1/kernel/Regularizer/mul/x:output:0)conv_1/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: v
-conv_1/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv_1_10633*
_output_shapes
: *
dtype0�
conv_1/bias/Regularizer/L2LossL2Loss5conv_1/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/bias/Regularizer/mulMul&conv_1/bias/Regularizer/mul/x:output:0'conv_1/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: y
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp
fc_0_10637*
_output_shapes
:	�. *
dtype0�
fc_0/kernel/Regularizer/L2LossL2Loss5fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
fc_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/kernel/Regularizer/mulMul&fc_0/kernel/Regularizer/mul/x:output:0'fc_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: r
+fc_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp
fc_0_10639*
_output_shapes
: *
dtype0|
fc_0/bias/Regularizer/L2LossL2Loss3fc_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: `
fc_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/bias/Regularizer/mulMul$fc_0/bias/Regularizer/mul/x:output:0%fc_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: t
IdentityIdentity%fc_0/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:��������� �
NoOpNoOp^conv_0/StatefulPartitionedCall.^conv_0/bias/Regularizer/L2Loss/ReadVariableOp0^conv_0/kernel/Regularizer/L2Loss/ReadVariableOp^conv_1/StatefulPartitionedCall.^conv_1/bias/Regularizer/L2Loss/ReadVariableOp0^conv_1/kernel/Regularizer/L2Loss/ReadVariableOp^fc_0/StatefulPartitionedCall,^fc_0/bias/Regularizer/L2Loss/ReadVariableOp.^fc_0/kernel/Regularizer/L2Loss/ReadVariableOp'^gaussian_noise/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*:
_input_shapes)
':���������
': : : : : : 2@
conv_0/StatefulPartitionedCallconv_0/StatefulPartitionedCall2^
-conv_0/bias/Regularizer/L2Loss/ReadVariableOp-conv_0/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp2@
conv_1/StatefulPartitionedCallconv_1/StatefulPartitionedCall2^
-conv_1/bias/Regularizer/L2Loss/ReadVariableOp-conv_1/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp2<
fc_0/StatefulPartitionedCallfc_0/StatefulPartitionedCall2Z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOp+fc_0/bias/Regularizer/L2Loss/ReadVariableOp2^
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp2P
&gaussian_noise/StatefulPartitionedCall&gaussian_noise/StatefulPartitionedCall:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�	
�
,__inference_sequential_2_layer_call_fn_11485

inputs!
unknown: 
	unknown_0: #
	unknown_1:  
	unknown_2: 
	unknown_3:	�. 
	unknown_4: 
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *(
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_2_layer_call_and_return_conditional_losses_10667o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:��������� `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*:
_input_shapes)
':���������
': : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�	
�
0__inference_finetuning_model_layer_call_fn_11214

inputs!
unknown: 
	unknown_0: #
	unknown_1:  
	unknown_2: 
	unknown_3:	�. 
	unknown_4: 
	unknown_5: 

	unknown_6:

identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *T
fORM
K__inference_finetuning_model_layer_call_and_return_conditional_losses_10981o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������
`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:���������
': : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�
�
A__inference_conv_0_layer_call_and_return_conditional_losses_11736

inputs8
conv2d_readvariableop_resource: -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�Conv2D/ReadVariableOp�-conv_0/bias/Regularizer/L2Loss/ReadVariableOp�/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp|
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
: *
dtype0�
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������# *
paddingVALID*
strides
r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0}
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������# X
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:���������# �
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
: *
dtype0�
 conv_0/kernel/Regularizer/L2LossL2Loss7conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/kernel/Regularizer/mulMul(conv_0/kernel/Regularizer/mul/x:output:0)conv_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-conv_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_0/bias/Regularizer/L2LossL2Loss5conv_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/bias/Regularizer/mulMul&conv_0/bias/Regularizer/mul/x:output:0'conv_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: i
IdentityIdentityRelu:activations:0^NoOp*
T0*/
_output_shapes
:���������# �
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp.^conv_0/bias/Regularizer/L2Loss/ReadVariableOp0^conv_0/kernel/Regularizer/L2Loss/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������
': : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp2^
-conv_0/bias/Regularizer/L2Loss/ReadVariableOp-conv_0/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�	
�
0__inference_finetuning_model_layer_call_fn_11068
input_5!
unknown: 
	unknown_0: #
	unknown_1:  
	unknown_2: 
	unknown_3:	�. 
	unknown_4: 
	unknown_5: 

	unknown_6:

identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinput_5unknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
**
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *T
fORM
K__inference_finetuning_model_layer_call_and_return_conditional_losses_11049o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������
`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:���������
': : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:X T
/
_output_shapes
:���������
'
!
_user_specified_name	input_5
�
�
&__inference_conv_1_layer_call_fn_11745

inputs!
unknown:  
	unknown_0: 
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *J
fERC
A__inference_conv_1_layer_call_and_return_conditional_losses_10506w
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:��������� `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������# : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������# 
 
_user_specified_nameinputs
�
�
__inference_loss_fn_5_11857B
4fc_0_bias_regularizer_l2loss_readvariableop_resource: 
identity��+fc_0/bias/Regularizer/L2Loss/ReadVariableOp�
+fc_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp4fc_0_bias_regularizer_l2loss_readvariableop_resource*
_output_shapes
: *
dtype0|
fc_0/bias/Regularizer/L2LossL2Loss3fc_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: `
fc_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/bias/Regularizer/mulMul$fc_0/bias/Regularizer/mul/x:output:0%fc_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: [
IdentityIdentityfc_0/bias/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: t
NoOpNoOp,^fc_0/bias/Regularizer/L2Loss/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2Z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOp+fc_0/bias/Regularizer/L2Loss/ReadVariableOp
�
�
?__inference_fc_0_layer_call_and_return_conditional_losses_11803

inputs1
matmul_readvariableop_resource:	�. -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOp�+fc_0/bias/Regularizer/L2Loss/ReadVariableOp�-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpu
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�. *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:��������� �
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	�. *
dtype0�
fc_0/kernel/Regularizer/L2LossL2Loss5fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
fc_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/kernel/Regularizer/mulMul&fc_0/kernel/Regularizer/mul/x:output:0'fc_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
+fc_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0|
fc_0/bias/Regularizer/L2LossL2Loss3fc_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: `
fc_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/bias/Regularizer/mulMul$fc_0/bias/Regularizer/mul/x:output:0%fc_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:��������� �
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp,^fc_0/bias/Regularizer/L2Loss/ReadVariableOp.^fc_0/kernel/Regularizer/L2Loss/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:����������.: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2Z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOp+fc_0/bias/Regularizer/L2Loss/ReadVariableOp2^
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:P L
(
_output_shapes
:����������.
 
_user_specified_nameinputs
�	
�
,__inference_sequential_2_layer_call_fn_10682
gaussian_noise_input!
unknown: 
	unknown_0: #
	unknown_1:  
	unknown_2: 
	unknown_3:	�. 
	unknown_4: 
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallgaussian_noise_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *(
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_2_layer_call_and_return_conditional_losses_10667o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:��������� `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*:
_input_shapes)
':���������
': : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:e a
/
_output_shapes
:���������
'
.
_user_specified_namegaussian_noise_input
�
l
P__inference_random_color_affine_2_layer_call_and_return_conditional_losses_10400

images
identityV
IdentityIdentityimages*
T0*/
_output_shapes
:���������
'"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
':W S
/
_output_shapes
:���������
'
 
_user_specified_nameimages
�	
�
__inference_loss_fn_2_11830R
8conv_1_kernel_regularizer_l2loss_readvariableop_resource:  
identity��/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp�
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp8conv_1_kernel_regularizer_l2loss_readvariableop_resource*&
_output_shapes
:  *
dtype0�
 conv_1/kernel/Regularizer/L2LossL2Loss7conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/kernel/Regularizer/mulMul(conv_1/kernel/Regularizer/mul/x:output:0)conv_1/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: _
IdentityIdentity!conv_1/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: x
NoOpNoOp0^conv_1/kernel/Regularizer/L2Loss/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2b
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp
�
H
,__inference_sequential_3_layer_call_fn_11401

inputs
identity�
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_3_layer_call_and_return_conditional_losses_10421h
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:���������
'"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
':W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�	
�
,__inference_sequential_2_layer_call_fn_11502

inputs!
unknown: 
	unknown_0: #
	unknown_1:  
	unknown_2: 
	unknown_3:	�. 
	unknown_4: 
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *(
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_2_layer_call_and_return_conditional_losses_10729o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:��������� `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*:
_input_shapes)
':���������
': : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�H
�
G__inference_sequential_2_layer_call_and_return_conditional_losses_11560

inputs?
%conv_0_conv2d_readvariableop_resource: 4
&conv_0_biasadd_readvariableop_resource: ?
%conv_1_conv2d_readvariableop_resource:  4
&conv_1_biasadd_readvariableop_resource: 6
#fc_0_matmul_readvariableop_resource:	�. 2
$fc_0_biasadd_readvariableop_resource: 
identity��conv_0/BiasAdd/ReadVariableOp�conv_0/Conv2D/ReadVariableOp�-conv_0/bias/Regularizer/L2Loss/ReadVariableOp�/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp�conv_1/BiasAdd/ReadVariableOp�conv_1/Conv2D/ReadVariableOp�-conv_1/bias/Regularizer/L2Loss/ReadVariableOp�/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp�fc_0/BiasAdd/ReadVariableOp�fc_0/MatMul/ReadVariableOp�+fc_0/bias/Regularizer/L2Loss/ReadVariableOp�-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpX
gaussian_noise/ShapeShapeinputs*
T0*
_output_shapes
::��f
!gaussian_noise/random_normal/meanConst*
_output_shapes
: *
dtype0*
valueB
 *    h
#gaussian_noise/random_normal/stddevConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
1gaussian_noise/random_normal/RandomStandardNormalRandomStandardNormalgaussian_noise/Shape:output:0*
T0*/
_output_shapes
:���������
'*
dtype0�
 gaussian_noise/random_normal/mulMul:gaussian_noise/random_normal/RandomStandardNormal:output:0,gaussian_noise/random_normal/stddev:output:0*
T0*/
_output_shapes
:���������
'�
gaussian_noise/random_normalAddV2$gaussian_noise/random_normal/mul:z:0*gaussian_noise/random_normal/mean:output:0*
T0*/
_output_shapes
:���������
'
gaussian_noise/addAddV2inputs gaussian_noise/random_normal:z:0*
T0*/
_output_shapes
:���������
'�
conv_0/Conv2D/ReadVariableOpReadVariableOp%conv_0_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype0�
conv_0/Conv2DConv2Dgaussian_noise/add:z:0$conv_0/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������# *
paddingVALID*
strides
�
conv_0/BiasAdd/ReadVariableOpReadVariableOp&conv_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_0/BiasAddBiasAddconv_0/Conv2D:output:0%conv_0/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������# f
conv_0/ReluReluconv_0/BiasAdd:output:0*
T0*/
_output_shapes
:���������# �
conv_1/Conv2D/ReadVariableOpReadVariableOp%conv_1_conv2d_readvariableop_resource*&
_output_shapes
:  *
dtype0�
conv_1/Conv2DConv2Dconv_0/Relu:activations:0$conv_1/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:��������� *
paddingVALID*
strides
�
conv_1/BiasAdd/ReadVariableOpReadVariableOp&conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_1/BiasAddBiasAddconv_1/Conv2D:output:0%conv_1/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:��������� f
conv_1/ReluReluconv_1/BiasAdd:output:0*
T0*/
_output_shapes
:��������� ^
flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"����@  �
flatten/ReshapeReshapeconv_1/Relu:activations:0flatten/Const:output:0*
T0*(
_output_shapes
:����������.
fc_0/MatMul/ReadVariableOpReadVariableOp#fc_0_matmul_readvariableop_resource*
_output_shapes
:	�. *
dtype0�
fc_0/MatMulMatMulflatten/Reshape:output:0"fc_0/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� |
fc_0/BiasAdd/ReadVariableOpReadVariableOp$fc_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
fc_0/BiasAddBiasAddfc_0/MatMul:product:0#fc_0/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� Z
	fc_0/ReluRelufc_0/BiasAdd:output:0*
T0*'
_output_shapes
:��������� �
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp%conv_0_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype0�
 conv_0/kernel/Regularizer/L2LossL2Loss7conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/kernel/Regularizer/mulMul(conv_0/kernel/Regularizer/mul/x:output:0)conv_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-conv_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp&conv_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_0/bias/Regularizer/L2LossL2Loss5conv_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/bias/Regularizer/mulMul&conv_0/bias/Regularizer/mul/x:output:0'conv_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp%conv_1_conv2d_readvariableop_resource*&
_output_shapes
:  *
dtype0�
 conv_1/kernel/Regularizer/L2LossL2Loss7conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/kernel/Regularizer/mulMul(conv_1/kernel/Regularizer/mul/x:output:0)conv_1/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-conv_1/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp&conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_1/bias/Regularizer/L2LossL2Loss5conv_1/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/bias/Regularizer/mulMul&conv_1/bias/Regularizer/mul/x:output:0'conv_1/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp#fc_0_matmul_readvariableop_resource*
_output_shapes
:	�. *
dtype0�
fc_0/kernel/Regularizer/L2LossL2Loss5fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
fc_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/kernel/Regularizer/mulMul&fc_0/kernel/Regularizer/mul/x:output:0'fc_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
+fc_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp$fc_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0|
fc_0/bias/Regularizer/L2LossL2Loss3fc_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: `
fc_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/bias/Regularizer/mulMul$fc_0/bias/Regularizer/mul/x:output:0%fc_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: f
IdentityIdentityfc_0/Relu:activations:0^NoOp*
T0*'
_output_shapes
:��������� �
NoOpNoOp^conv_0/BiasAdd/ReadVariableOp^conv_0/Conv2D/ReadVariableOp.^conv_0/bias/Regularizer/L2Loss/ReadVariableOp0^conv_0/kernel/Regularizer/L2Loss/ReadVariableOp^conv_1/BiasAdd/ReadVariableOp^conv_1/Conv2D/ReadVariableOp.^conv_1/bias/Regularizer/L2Loss/ReadVariableOp0^conv_1/kernel/Regularizer/L2Loss/ReadVariableOp^fc_0/BiasAdd/ReadVariableOp^fc_0/MatMul/ReadVariableOp,^fc_0/bias/Regularizer/L2Loss/ReadVariableOp.^fc_0/kernel/Regularizer/L2Loss/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*:
_input_shapes)
':���������
': : : : : : 2>
conv_0/BiasAdd/ReadVariableOpconv_0/BiasAdd/ReadVariableOp2<
conv_0/Conv2D/ReadVariableOpconv_0/Conv2D/ReadVariableOp2^
-conv_0/bias/Regularizer/L2Loss/ReadVariableOp-conv_0/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp2>
conv_1/BiasAdd/ReadVariableOpconv_1/BiasAdd/ReadVariableOp2<
conv_1/Conv2D/ReadVariableOpconv_1/Conv2D/ReadVariableOp2^
-conv_1/bias/Regularizer/L2Loss/ReadVariableOp-conv_1/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp2:
fc_0/BiasAdd/ReadVariableOpfc_0/BiasAdd/ReadVariableOp28
fc_0/MatMul/ReadVariableOpfc_0/MatMul/ReadVariableOp2Z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOp+fc_0/bias/Regularizer/L2Loss/ReadVariableOp2^
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�3
�
K__inference_finetuning_model_layer_call_and_return_conditional_losses_11049

inputs,
sequential_2_11006:  
sequential_2_11008: ,
sequential_2_11010:   
sequential_2_11012: %
sequential_2_11014:	�.  
sequential_2_11016: 
dense_3_11019: 

dense_3_11021:

identity��-conv_0/bias/Regularizer/L2Loss/ReadVariableOp�/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp�-conv_1/bias/Regularizer/L2Loss/ReadVariableOp�/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp�dense_3/StatefulPartitionedCall�+fc_0/bias/Regularizer/L2Loss/ReadVariableOp�-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp�$sequential_2/StatefulPartitionedCall�
sequential_3/PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_3_layer_call_and_return_conditional_losses_10421�
$sequential_2/StatefulPartitionedCallStatefulPartitionedCall%sequential_3/PartitionedCall:output:0sequential_2_11006sequential_2_11008sequential_2_11010sequential_2_11012sequential_2_11014sequential_2_11016*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *(
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_2_layer_call_and_return_conditional_losses_10729�
dense_3/StatefulPartitionedCallStatefulPartitionedCall-sequential_2/StatefulPartitionedCall:output:0dense_3_11019dense_3_11021*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_3_layer_call_and_return_conditional_losses_10853�
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_11006*&
_output_shapes
: *
dtype0�
 conv_0/kernel/Regularizer/L2LossL2Loss7conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/kernel/Regularizer/mulMul(conv_0/kernel/Regularizer/mul/x:output:0)conv_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: |
-conv_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_11008*
_output_shapes
: *
dtype0�
conv_0/bias/Regularizer/L2LossL2Loss5conv_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/bias/Regularizer/mulMul&conv_0/bias/Regularizer/mul/x:output:0'conv_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_11010*&
_output_shapes
:  *
dtype0�
 conv_1/kernel/Regularizer/L2LossL2Loss7conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/kernel/Regularizer/mulMul(conv_1/kernel/Regularizer/mul/x:output:0)conv_1/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: |
-conv_1/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_11012*
_output_shapes
: *
dtype0�
conv_1/bias/Regularizer/L2LossL2Loss5conv_1/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/bias/Regularizer/mulMul&conv_1/bias/Regularizer/mul/x:output:0'conv_1/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_11014*
_output_shapes
:	�. *
dtype0�
fc_0/kernel/Regularizer/L2LossL2Loss5fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
fc_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/kernel/Regularizer/mulMul&fc_0/kernel/Regularizer/mul/x:output:0'fc_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpsequential_2_11016*
_output_shapes
: *
dtype0|
fc_0/bias/Regularizer/L2LossL2Loss3fc_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: `
fc_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/bias/Regularizer/mulMul$fc_0/bias/Regularizer/mul/x:output:0%fc_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: w
IdentityIdentity(dense_3/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������
�
NoOpNoOp.^conv_0/bias/Regularizer/L2Loss/ReadVariableOp0^conv_0/kernel/Regularizer/L2Loss/ReadVariableOp.^conv_1/bias/Regularizer/L2Loss/ReadVariableOp0^conv_1/kernel/Regularizer/L2Loss/ReadVariableOp ^dense_3/StatefulPartitionedCall,^fc_0/bias/Regularizer/L2Loss/ReadVariableOp.^fc_0/kernel/Regularizer/L2Loss/ReadVariableOp%^sequential_2/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:���������
': : : : : : : : 2^
-conv_0/bias/Regularizer/L2Loss/ReadVariableOp-conv_0/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp2^
-conv_1/bias/Regularizer/L2Loss/ReadVariableOp-conv_1/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2Z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOp+fc_0/bias/Regularizer/L2Loss/ReadVariableOp2^
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp2L
$sequential_2/StatefulPartitionedCall$sequential_2/StatefulPartitionedCall:W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�
�
G__inference_sequential_3_layer_call_and_return_conditional_losses_10394
input_6
identity��-random_color_affine_2/StatefulPartitionedCall�
-random_color_affine_2/StatefulPartitionedCallStatefulPartitionedCallinput_6*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *Y
fTRR
P__inference_random_color_affine_2_layer_call_and_return_conditional_losses_10391�
IdentityIdentity6random_color_affine_2/StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:���������
'v
NoOpNoOp.^random_color_affine_2/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
'2^
-random_color_affine_2/StatefulPartitionedCall-random_color_affine_2/StatefulPartitionedCall:X T
/
_output_shapes
:���������
'
!
_user_specified_name	input_6
�
�
A__inference_conv_1_layer_call_and_return_conditional_losses_11764

inputs8
conv2d_readvariableop_resource:  -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�Conv2D/ReadVariableOp�-conv_1/bias/Regularizer/L2Loss/ReadVariableOp�/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp|
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:  *
dtype0�
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:��������� *
paddingVALID*
strides
r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0}
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:��������� X
ReluReluBiasAdd:output:0*
T0*/
_output_shapes
:��������� �
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:  *
dtype0�
 conv_1/kernel/Regularizer/L2LossL2Loss7conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/kernel/Regularizer/mulMul(conv_1/kernel/Regularizer/mul/x:output:0)conv_1/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
-conv_1/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_1/bias/Regularizer/L2LossL2Loss5conv_1/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/bias/Regularizer/mulMul&conv_1/bias/Regularizer/mul/x:output:0'conv_1/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: i
IdentityIdentityRelu:activations:0^NoOp*
T0*/
_output_shapes
:��������� �
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp.^conv_1/bias/Regularizer/L2Loss/ReadVariableOp0^conv_1/kernel/Regularizer/L2Loss/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:���������# : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp2^
-conv_1/bias/Regularizer/L2Loss/ReadVariableOp-conv_1/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:W S
/
_output_shapes
:���������# 
 
_user_specified_nameinputs
�
d
G__inference_sequential_3_layer_call_and_return_conditional_losses_10403
input_6
identity�
%random_color_affine_2/PartitionedCallPartitionedCallinput_6*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *Y
fTRR
P__inference_random_color_affine_2_layer_call_and_return_conditional_losses_10400~
IdentityIdentity.random_color_affine_2/PartitionedCall:output:0*
T0*/
_output_shapes
:���������
'"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
':X T
/
_output_shapes
:���������
'
!
_user_specified_name	input_6
�
l
P__inference_random_color_affine_2_layer_call_and_return_conditional_losses_11683

images
identityV
IdentityIdentityimages*
T0*/
_output_shapes
:���������
'"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
':W S
/
_output_shapes
:���������
'
 
_user_specified_nameimages
�8
�
G__inference_sequential_2_layer_call_and_return_conditional_losses_10619
gaussian_noise_input&
conv_0_10578: 
conv_0_10580: &
conv_1_10583:  
conv_1_10585: 

fc_0_10589:	�. 

fc_0_10591: 
identity��conv_0/StatefulPartitionedCall�-conv_0/bias/Regularizer/L2Loss/ReadVariableOp�/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp�conv_1/StatefulPartitionedCall�-conv_1/bias/Regularizer/L2Loss/ReadVariableOp�/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp�fc_0/StatefulPartitionedCall�+fc_0/bias/Regularizer/L2Loss/ReadVariableOp�-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp�
gaussian_noise/PartitionedCallPartitionedCallgaussian_noise_input*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������
'* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *R
fMRK
I__inference_gaussian_noise_layer_call_and_return_conditional_losses_10576�
conv_0/StatefulPartitionedCallStatefulPartitionedCall'gaussian_noise/PartitionedCall:output:0conv_0_10578conv_0_10580*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:���������# *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *J
fERC
A__inference_conv_0_layer_call_and_return_conditional_losses_10481�
conv_1/StatefulPartitionedCallStatefulPartitionedCall'conv_0/StatefulPartitionedCall:output:0conv_1_10583conv_1_10585*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *J
fERC
A__inference_conv_1_layer_call_and_return_conditional_losses_10506�
flatten/PartitionedCallPartitionedCall'conv_1/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:����������.* 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_flatten_layer_call_and_return_conditional_losses_10518�
fc_0/StatefulPartitionedCallStatefulPartitionedCall flatten/PartitionedCall:output:0
fc_0_10589
fc_0_10591*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_fc_0_layer_call_and_return_conditional_losses_10539�
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv_0_10578*&
_output_shapes
: *
dtype0�
 conv_0/kernel/Regularizer/L2LossL2Loss7conv_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/kernel/Regularizer/mulMul(conv_0/kernel/Regularizer/mul/x:output:0)conv_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: v
-conv_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv_0_10580*
_output_shapes
: *
dtype0�
conv_0/bias/Regularizer/L2LossL2Loss5conv_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/bias/Regularizer/mulMul&conv_0/bias/Regularizer/mul/x:output:0'conv_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: �
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv_1_10583*&
_output_shapes
:  *
dtype0�
 conv_1/kernel/Regularizer/L2LossL2Loss7conv_1/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: d
conv_1/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/kernel/Regularizer/mulMul(conv_1/kernel/Regularizer/mul/x:output:0)conv_1/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: v
-conv_1/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOpconv_1_10585*
_output_shapes
: *
dtype0�
conv_1/bias/Regularizer/L2LossL2Loss5conv_1/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_1/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_1/bias/Regularizer/mulMul&conv_1/bias/Regularizer/mul/x:output:0'conv_1/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: y
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOpReadVariableOp
fc_0_10589*
_output_shapes
:	�. *
dtype0�
fc_0/kernel/Regularizer/L2LossL2Loss5fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
fc_0/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/kernel/Regularizer/mulMul&fc_0/kernel/Regularizer/mul/x:output:0'fc_0/kernel/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: r
+fc_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp
fc_0_10591*
_output_shapes
: *
dtype0|
fc_0/bias/Regularizer/L2LossL2Loss3fc_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: `
fc_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
fc_0/bias/Regularizer/mulMul$fc_0/bias/Regularizer/mul/x:output:0%fc_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: t
IdentityIdentity%fc_0/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:��������� �
NoOpNoOp^conv_0/StatefulPartitionedCall.^conv_0/bias/Regularizer/L2Loss/ReadVariableOp0^conv_0/kernel/Regularizer/L2Loss/ReadVariableOp^conv_1/StatefulPartitionedCall.^conv_1/bias/Regularizer/L2Loss/ReadVariableOp0^conv_1/kernel/Regularizer/L2Loss/ReadVariableOp^fc_0/StatefulPartitionedCall,^fc_0/bias/Regularizer/L2Loss/ReadVariableOp.^fc_0/kernel/Regularizer/L2Loss/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*:
_input_shapes)
':���������
': : : : : : 2@
conv_0/StatefulPartitionedCallconv_0/StatefulPartitionedCall2^
-conv_0/bias/Regularizer/L2Loss/ReadVariableOp-conv_0/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp/conv_0/kernel/Regularizer/L2Loss/ReadVariableOp2@
conv_1/StatefulPartitionedCallconv_1/StatefulPartitionedCall2^
-conv_1/bias/Regularizer/L2Loss/ReadVariableOp-conv_1/bias/Regularizer/L2Loss/ReadVariableOp2b
/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp/conv_1/kernel/Regularizer/L2Loss/ReadVariableOp2<
fc_0/StatefulPartitionedCallfc_0/StatefulPartitionedCall2Z
+fc_0/bias/Regularizer/L2Loss/ReadVariableOp+fc_0/bias/Regularizer/L2Loss/ReadVariableOp2^
-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp-fc_0/kernel/Regularizer/L2Loss/ReadVariableOp:e a
/
_output_shapes
:���������
'
.
_user_specified_namegaussian_noise_input
�
�
'__inference_dense_3_layer_call_fn_11620

inputs
unknown: 

	unknown_0:

identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������
*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_dense_3_layer_call_and_return_conditional_losses_10853o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������
`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:��������� 
 
_user_specified_nameinputs
�
�
__inference_loss_fn_1_11821D
6conv_0_bias_regularizer_l2loss_readvariableop_resource: 
identity��-conv_0/bias/Regularizer/L2Loss/ReadVariableOp�
-conv_0/bias/Regularizer/L2Loss/ReadVariableOpReadVariableOp6conv_0_bias_regularizer_l2loss_readvariableop_resource*
_output_shapes
: *
dtype0�
conv_0/bias/Regularizer/L2LossL2Loss5conv_0/bias/Regularizer/L2Loss/ReadVariableOp:value:0*
T0*
_output_shapes
: b
conv_0/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *o;�
conv_0/bias/Regularizer/mulMul&conv_0/bias/Regularizer/mul/x:output:0'conv_0/bias/Regularizer/L2Loss:output:0*
T0*
_output_shapes
: ]
IdentityIdentityconv_0/bias/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: v
NoOpNoOp.^conv_0/bias/Regularizer/L2Loss/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2^
-conv_0/bias/Regularizer/L2Loss/ReadVariableOp-conv_0/bias/Regularizer/L2Loss/ReadVariableOp
�8
�	
 __inference__wrapped_model_10348
input_5]
Cfinetuning_model_sequential_2_conv_0_conv2d_readvariableop_resource: R
Dfinetuning_model_sequential_2_conv_0_biasadd_readvariableop_resource: ]
Cfinetuning_model_sequential_2_conv_1_conv2d_readvariableop_resource:  R
Dfinetuning_model_sequential_2_conv_1_biasadd_readvariableop_resource: T
Afinetuning_model_sequential_2_fc_0_matmul_readvariableop_resource:	�. P
Bfinetuning_model_sequential_2_fc_0_biasadd_readvariableop_resource: I
7finetuning_model_dense_3_matmul_readvariableop_resource: 
F
8finetuning_model_dense_3_biasadd_readvariableop_resource:

identity��/finetuning_model/dense_3/BiasAdd/ReadVariableOp�.finetuning_model/dense_3/MatMul/ReadVariableOp�;finetuning_model/sequential_2/conv_0/BiasAdd/ReadVariableOp�:finetuning_model/sequential_2/conv_0/Conv2D/ReadVariableOp�;finetuning_model/sequential_2/conv_1/BiasAdd/ReadVariableOp�:finetuning_model/sequential_2/conv_1/Conv2D/ReadVariableOp�9finetuning_model/sequential_2/fc_0/BiasAdd/ReadVariableOp�8finetuning_model/sequential_2/fc_0/MatMul/ReadVariableOp�
:finetuning_model/sequential_2/conv_0/Conv2D/ReadVariableOpReadVariableOpCfinetuning_model_sequential_2_conv_0_conv2d_readvariableop_resource*&
_output_shapes
: *
dtype0�
+finetuning_model/sequential_2/conv_0/Conv2DConv2Dinput_5Bfinetuning_model/sequential_2/conv_0/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������# *
paddingVALID*
strides
�
;finetuning_model/sequential_2/conv_0/BiasAdd/ReadVariableOpReadVariableOpDfinetuning_model_sequential_2_conv_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
,finetuning_model/sequential_2/conv_0/BiasAddBiasAdd4finetuning_model/sequential_2/conv_0/Conv2D:output:0Cfinetuning_model/sequential_2/conv_0/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:���������# �
)finetuning_model/sequential_2/conv_0/ReluRelu5finetuning_model/sequential_2/conv_0/BiasAdd:output:0*
T0*/
_output_shapes
:���������# �
:finetuning_model/sequential_2/conv_1/Conv2D/ReadVariableOpReadVariableOpCfinetuning_model_sequential_2_conv_1_conv2d_readvariableop_resource*&
_output_shapes
:  *
dtype0�
+finetuning_model/sequential_2/conv_1/Conv2DConv2D7finetuning_model/sequential_2/conv_0/Relu:activations:0Bfinetuning_model/sequential_2/conv_1/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:��������� *
paddingVALID*
strides
�
;finetuning_model/sequential_2/conv_1/BiasAdd/ReadVariableOpReadVariableOpDfinetuning_model_sequential_2_conv_1_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
,finetuning_model/sequential_2/conv_1/BiasAddBiasAdd4finetuning_model/sequential_2/conv_1/Conv2D:output:0Cfinetuning_model/sequential_2/conv_1/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:��������� �
)finetuning_model/sequential_2/conv_1/ReluRelu5finetuning_model/sequential_2/conv_1/BiasAdd:output:0*
T0*/
_output_shapes
:��������� |
+finetuning_model/sequential_2/flatten/ConstConst*
_output_shapes
:*
dtype0*
valueB"����@  �
-finetuning_model/sequential_2/flatten/ReshapeReshape7finetuning_model/sequential_2/conv_1/Relu:activations:04finetuning_model/sequential_2/flatten/Const:output:0*
T0*(
_output_shapes
:����������.�
8finetuning_model/sequential_2/fc_0/MatMul/ReadVariableOpReadVariableOpAfinetuning_model_sequential_2_fc_0_matmul_readvariableop_resource*
_output_shapes
:	�. *
dtype0�
)finetuning_model/sequential_2/fc_0/MatMulMatMul6finetuning_model/sequential_2/flatten/Reshape:output:0@finetuning_model/sequential_2/fc_0/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� �
9finetuning_model/sequential_2/fc_0/BiasAdd/ReadVariableOpReadVariableOpBfinetuning_model_sequential_2_fc_0_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
*finetuning_model/sequential_2/fc_0/BiasAddBiasAdd3finetuning_model/sequential_2/fc_0/MatMul:product:0Afinetuning_model/sequential_2/fc_0/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� �
'finetuning_model/sequential_2/fc_0/ReluRelu3finetuning_model/sequential_2/fc_0/BiasAdd:output:0*
T0*'
_output_shapes
:��������� �
.finetuning_model/dense_3/MatMul/ReadVariableOpReadVariableOp7finetuning_model_dense_3_matmul_readvariableop_resource*
_output_shapes

: 
*
dtype0�
finetuning_model/dense_3/MatMulMatMul5finetuning_model/sequential_2/fc_0/Relu:activations:06finetuning_model/dense_3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
�
/finetuning_model/dense_3/BiasAdd/ReadVariableOpReadVariableOp8finetuning_model_dense_3_biasadd_readvariableop_resource*
_output_shapes
:
*
dtype0�
 finetuning_model/dense_3/BiasAddBiasAdd)finetuning_model/dense_3/MatMul:product:07finetuning_model/dense_3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������
x
IdentityIdentity)finetuning_model/dense_3/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:���������
�
NoOpNoOp0^finetuning_model/dense_3/BiasAdd/ReadVariableOp/^finetuning_model/dense_3/MatMul/ReadVariableOp<^finetuning_model/sequential_2/conv_0/BiasAdd/ReadVariableOp;^finetuning_model/sequential_2/conv_0/Conv2D/ReadVariableOp<^finetuning_model/sequential_2/conv_1/BiasAdd/ReadVariableOp;^finetuning_model/sequential_2/conv_1/Conv2D/ReadVariableOp:^finetuning_model/sequential_2/fc_0/BiasAdd/ReadVariableOp9^finetuning_model/sequential_2/fc_0/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:���������
': : : : : : : : 2b
/finetuning_model/dense_3/BiasAdd/ReadVariableOp/finetuning_model/dense_3/BiasAdd/ReadVariableOp2`
.finetuning_model/dense_3/MatMul/ReadVariableOp.finetuning_model/dense_3/MatMul/ReadVariableOp2z
;finetuning_model/sequential_2/conv_0/BiasAdd/ReadVariableOp;finetuning_model/sequential_2/conv_0/BiasAdd/ReadVariableOp2x
:finetuning_model/sequential_2/conv_0/Conv2D/ReadVariableOp:finetuning_model/sequential_2/conv_0/Conv2D/ReadVariableOp2z
;finetuning_model/sequential_2/conv_1/BiasAdd/ReadVariableOp;finetuning_model/sequential_2/conv_1/BiasAdd/ReadVariableOp2x
:finetuning_model/sequential_2/conv_1/Conv2D/ReadVariableOp:finetuning_model/sequential_2/conv_1/Conv2D/ReadVariableOp2v
9finetuning_model/sequential_2/fc_0/BiasAdd/ReadVariableOp9finetuning_model/sequential_2/fc_0/BiasAdd/ReadVariableOp2t
8finetuning_model/sequential_2/fc_0/MatMul/ReadVariableOp8finetuning_model/sequential_2/fc_0/MatMul/ReadVariableOp:X T
/
_output_shapes
:���������
'
!
_user_specified_name	input_5
�	
�
,__inference_sequential_2_layer_call_fn_10744
gaussian_noise_input!
unknown: 
	unknown_0: #
	unknown_1:  
	unknown_2: 
	unknown_3:	�. 
	unknown_4: 
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallgaussian_noise_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4*
Tin
	2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *(
_read_only_resource_inputs

*-
config_proto

CPU

GPU 2J 8� *P
fKRI
G__inference_sequential_2_layer_call_and_return_conditional_losses_10729o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:��������� `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime*:
_input_shapes)
':���������
': : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:e a
/
_output_shapes
:���������
'
.
_user_specified_namegaussian_noise_input
�	
h
I__inference_gaussian_noise_layer_call_and_return_conditional_losses_11704

inputs
identity�I
ShapeShapeinputs*
T0*
_output_shapes
::��W
random_normal/meanConst*
_output_shapes
: *
dtype0*
valueB
 *    Y
random_normal/stddevConst*
_output_shapes
: *
dtype0*
valueB
 *
�#<�
"random_normal/RandomStandardNormalRandomStandardNormalShape:output:0*
T0*/
_output_shapes
:���������
'*
dtype0�
random_normal/mulMul+random_normal/RandomStandardNormal:output:0random_normal/stddev:output:0*
T0*/
_output_shapes
:���������
'�
random_normalAddV2random_normal/mul:z:0random_normal/mean:output:0*
T0*/
_output_shapes
:���������
'a
addAddV2inputsrandom_normal:z:0*
T0*/
_output_shapes
:���������
'W
IdentityIdentityadd:z:0*
T0*/
_output_shapes
:���������
'"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
':W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs
�
c
G__inference_sequential_3_layer_call_and_return_conditional_losses_11444

inputs
identityV
IdentityIdentityinputs*
T0*/
_output_shapes
:���������
'"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:���������
':W S
/
_output_shapes
:���������
'
 
_user_specified_nameinputs"�
L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*�
serving_default�
C
input_58
serving_default_input_5:0���������
';
dense_30
StatefulPartitionedCall:0���������
tensorflow/serving/predict:��
�
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*	&call_and_return_all_conditional_losses

_default_save_signature
	optimizer

signatures"
_tf_keras_sequential
�
layer-0
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses"
_tf_keras_sequential
�
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer-3
layer_with_weights-2
layer-4
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses"
_tf_keras_sequential
�
	variables
 trainable_variables
!regularization_losses
"	keras_api
#__call__
*$&call_and_return_all_conditional_losses

%kernel
&bias"
_tf_keras_layer
X
'0
(1
)2
*3
+4
,5
%6
&7"
trackable_list_wrapper
X
'0
(1
)2
*3
+4
,5
%6
&7"
trackable_list_wrapper
 "
trackable_list_wrapper
�
-non_trainable_variables

.layers
/metrics
0layer_regularization_losses
1layer_metrics
	variables
trainable_variables
regularization_losses
__call__

_default_save_signature
*	&call_and_return_all_conditional_losses
&	"call_and_return_conditional_losses"
_generic_user_object
�
2trace_0
3trace_1
4trace_2
5trace_32�
0__inference_finetuning_model_layer_call_fn_11000
0__inference_finetuning_model_layer_call_fn_11068
0__inference_finetuning_model_layer_call_fn_11214
0__inference_finetuning_model_layer_call_fn_11235�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z2trace_0z3trace_1z4trace_2z5trace_3
�
6trace_0
7trace_1
8trace_2
9trace_32�
K__inference_finetuning_model_layer_call_and_return_conditional_losses_10884
K__inference_finetuning_model_layer_call_and_return_conditional_losses_10931
K__inference_finetuning_model_layer_call_and_return_conditional_losses_11334
K__inference_finetuning_model_layer_call_and_return_conditional_losses_11391�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z6trace_0z7trace_1z8trace_2z9trace_3
�B�
 __inference__wrapped_model_10348input_5"�
���
FullArgSpec
args� 
varargsjargs
varkwjkwargs
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�
:
_variables
;_iterations
<_learning_rate
=_index_dict
>
_momentums
?_velocities
@_update_step_xla"
experimentalOptimizer
,
Aserving_default"
signature_map
�
B	variables
Ctrainable_variables
Dregularization_losses
E	keras_api
F__call__
*G&call_and_return_all_conditional_losses"
_tf_keras_layer
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
Hnon_trainable_variables

Ilayers
Jmetrics
Klayer_regularization_losses
Llayer_metrics
	variables
trainable_variables
regularization_losses
__call__
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses"
_generic_user_object
�
Mtrace_0
Ntrace_1
Otrace_2
Ptrace_32�
,__inference_sequential_3_layer_call_fn_10414
,__inference_sequential_3_layer_call_fn_10424
,__inference_sequential_3_layer_call_fn_11396
,__inference_sequential_3_layer_call_fn_11401�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 zMtrace_0zNtrace_1zOtrace_2zPtrace_3
�
Qtrace_0
Rtrace_1
Strace_2
Ttrace_32�
G__inference_sequential_3_layer_call_and_return_conditional_losses_10394
G__inference_sequential_3_layer_call_and_return_conditional_losses_10403
G__inference_sequential_3_layer_call_and_return_conditional_losses_11440
G__inference_sequential_3_layer_call_and_return_conditional_losses_11444�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 zQtrace_0zRtrace_1zStrace_2zTtrace_3
�
U	variables
Vtrainable_variables
Wregularization_losses
X	keras_api
Y__call__
*Z&call_and_return_all_conditional_losses
[_random_generator"
_tf_keras_layer
�
\	variables
]trainable_variables
^regularization_losses
_	keras_api
`__call__
*a&call_and_return_all_conditional_losses

'kernel
(bias
 b_jit_compiled_convolution_op"
_tf_keras_layer
�
c	variables
dtrainable_variables
eregularization_losses
f	keras_api
g__call__
*h&call_and_return_all_conditional_losses

)kernel
*bias
 i_jit_compiled_convolution_op"
_tf_keras_layer
�
j	variables
ktrainable_variables
lregularization_losses
m	keras_api
n__call__
*o&call_and_return_all_conditional_losses"
_tf_keras_layer
�
p	variables
qtrainable_variables
rregularization_losses
s	keras_api
t__call__
*u&call_and_return_all_conditional_losses

+kernel
,bias"
_tf_keras_layer
J
'0
(1
)2
*3
+4
,5"
trackable_list_wrapper
J
'0
(1
)2
*3
+4
,5"
trackable_list_wrapper
J
v0
w1
x2
y3
z4
{5"
trackable_list_wrapper
�
|non_trainable_variables

}layers
~metrics
layer_regularization_losses
�layer_metrics
	variables
trainable_variables
regularization_losses
__call__
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses"
_generic_user_object
�
�trace_0
�trace_1
�trace_2
�trace_32�
,__inference_sequential_2_layer_call_fn_10682
,__inference_sequential_2_layer_call_fn_10744
,__inference_sequential_2_layer_call_fn_11485
,__inference_sequential_2_layer_call_fn_11502�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0z�trace_1z�trace_2z�trace_3
�
�trace_0
�trace_1
�trace_2
�trace_32�
G__inference_sequential_2_layer_call_and_return_conditional_losses_10570
G__inference_sequential_2_layer_call_and_return_conditional_losses_10619
G__inference_sequential_2_layer_call_and_return_conditional_losses_11560
G__inference_sequential_2_layer_call_and_return_conditional_losses_11611�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0z�trace_1z�trace_2z�trace_3
.
%0
&1"
trackable_list_wrapper
.
%0
&1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
	variables
 trainable_variables
!regularization_losses
#__call__
*$&call_and_return_all_conditional_losses
&$"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
'__inference_dense_3_layer_call_fn_11620�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
B__inference_dense_3_layer_call_and_return_conditional_losses_11630�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
 : 
2dense_3/kernel
:
2dense_3/bias
':% 2conv_0/kernel
: 2conv_0/bias
':%  2conv_1/kernel
: 2conv_1/bias
:	�. 2fc_0/kernel
: 2	fc_0/bias
 "
trackable_list_wrapper
5
0
1
2"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
0__inference_finetuning_model_layer_call_fn_11000input_5"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
0__inference_finetuning_model_layer_call_fn_11068input_5"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
0__inference_finetuning_model_layer_call_fn_11214inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
0__inference_finetuning_model_layer_call_fn_11235inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
K__inference_finetuning_model_layer_call_and_return_conditional_losses_10884input_5"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
K__inference_finetuning_model_layer_call_and_return_conditional_losses_10931input_5"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
K__inference_finetuning_model_layer_call_and_return_conditional_losses_11334inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
K__inference_finetuning_model_layer_call_and_return_conditional_losses_11391inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�
;0
�1
�2
�3
�4
�5
�6
�7
�8
�9
�10
�11
�12
�13
�14
�15
�16"
trackable_list_wrapper
:	 2	iteration
: 2learning_rate
 "
trackable_dict_wrapper
`
�0
�1
�2
�3
�4
�5
�6
�7"
trackable_list_wrapper
`
�0
�1
�2
�3
�4
�5
�6
�7"
trackable_list_wrapper
�2��
���
FullArgSpec*
args"�

jgradient

jvariable
jkey
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 0
�B�
#__inference_signature_wrapper_11169input_5"�
���
FullArgSpec
args� 
varargs
 
varkwjkwargs
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
B	variables
Ctrainable_variables
Dregularization_losses
F__call__
*G&call_and_return_all_conditional_losses
&G"call_and_return_conditional_losses"
_generic_user_object
�
�trace_0
�trace_12�
5__inference_random_color_affine_2_layer_call_fn_11635
5__inference_random_color_affine_2_layer_call_fn_11640�
���
FullArgSpec!
args�
jimages

jtraining
varargs
 
varkw
 
defaults�
p

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0z�trace_1
�
�trace_0
�trace_12�
P__inference_random_color_affine_2_layer_call_and_return_conditional_losses_11679
P__inference_random_color_affine_2_layer_call_and_return_conditional_losses_11683�
���
FullArgSpec!
args�
jimages

jtraining
varargs
 
varkw
 
defaults�
p

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0z�trace_1
 "
trackable_list_wrapper
'
0"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
,__inference_sequential_3_layer_call_fn_10414input_6"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
,__inference_sequential_3_layer_call_fn_10424input_6"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
,__inference_sequential_3_layer_call_fn_11396inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
,__inference_sequential_3_layer_call_fn_11401inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
G__inference_sequential_3_layer_call_and_return_conditional_losses_10394input_6"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
G__inference_sequential_3_layer_call_and_return_conditional_losses_10403input_6"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
G__inference_sequential_3_layer_call_and_return_conditional_losses_11440inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
G__inference_sequential_3_layer_call_and_return_conditional_losses_11444inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
U	variables
Vtrainable_variables
Wregularization_losses
Y__call__
*Z&call_and_return_all_conditional_losses
&Z"call_and_return_conditional_losses"
_generic_user_object
�
�trace_0
�trace_12�
.__inference_gaussian_noise_layer_call_fn_11688
.__inference_gaussian_noise_layer_call_fn_11693�
���
FullArgSpec!
args�
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0z�trace_1
�
�trace_0
�trace_12�
I__inference_gaussian_noise_layer_call_and_return_conditional_losses_11704
I__inference_gaussian_noise_layer_call_and_return_conditional_losses_11708�
���
FullArgSpec!
args�
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0z�trace_1
"
_generic_user_object
.
'0
(1"
trackable_list_wrapper
.
'0
(1"
trackable_list_wrapper
.
v0
w1"
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
\	variables
]trainable_variables
^regularization_losses
`__call__
*a&call_and_return_all_conditional_losses
&a"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
&__inference_conv_0_layer_call_fn_11717�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
A__inference_conv_0_layer_call_and_return_conditional_losses_11736�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�2��
���
FullArgSpec
args�
jinputs
jkernel
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 0
.
)0
*1"
trackable_list_wrapper
.
)0
*1"
trackable_list_wrapper
.
x0
y1"
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
c	variables
dtrainable_variables
eregularization_losses
g__call__
*h&call_and_return_all_conditional_losses
&h"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
&__inference_conv_1_layer_call_fn_11745�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
A__inference_conv_1_layer_call_and_return_conditional_losses_11764�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�2��
���
FullArgSpec
args�
jinputs
jkernel
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 0
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
j	variables
ktrainable_variables
lregularization_losses
n__call__
*o&call_and_return_all_conditional_losses
&o"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
'__inference_flatten_layer_call_fn_11769�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
B__inference_flatten_layer_call_and_return_conditional_losses_11775�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
.
+0
,1"
trackable_list_wrapper
.
+0
,1"
trackable_list_wrapper
.
z0
{1"
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
p	variables
qtrainable_variables
rregularization_losses
t__call__
*u&call_and_return_all_conditional_losses
&u"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
$__inference_fc_0_layer_call_fn_11784�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
?__inference_fc_0_layer_call_and_return_conditional_losses_11803�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
__inference_loss_fn_0_11812�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� z�trace_0
�
�trace_02�
__inference_loss_fn_1_11821�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� z�trace_0
�
�trace_02�
__inference_loss_fn_2_11830�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� z�trace_0
�
�trace_02�
__inference_loss_fn_3_11839�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� z�trace_0
�
�trace_02�
__inference_loss_fn_4_11848�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� z�trace_0
�
�trace_02�
__inference_loss_fn_5_11857�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� z�trace_0
 "
trackable_list_wrapper
C
0
1
2
3
4"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
,__inference_sequential_2_layer_call_fn_10682gaussian_noise_input"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
,__inference_sequential_2_layer_call_fn_10744gaussian_noise_input"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
,__inference_sequential_2_layer_call_fn_11485inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
,__inference_sequential_2_layer_call_fn_11502inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
G__inference_sequential_2_layer_call_and_return_conditional_losses_10570gaussian_noise_input"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
G__inference_sequential_2_layer_call_and_return_conditional_losses_10619gaussian_noise_input"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
G__inference_sequential_2_layer_call_and_return_conditional_losses_11560inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
G__inference_sequential_2_layer_call_and_return_conditional_losses_11611inputs"�
���
FullArgSpec)
args!�
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
'__inference_dense_3_layer_call_fn_11620inputs"�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
B__inference_dense_3_layer_call_and_return_conditional_losses_11630inputs"�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
R
�	variables
�	keras_api

�total

�count"
_tf_keras_metric
c
�	variables
�	keras_api

�total

�count
�
_fn_kwargs"
_tf_keras_metric
,:* 2Adam/m/conv_0/kernel
,:* 2Adam/v/conv_0/kernel
: 2Adam/m/conv_0/bias
: 2Adam/v/conv_0/bias
,:*  2Adam/m/conv_1/kernel
,:*  2Adam/v/conv_1/kernel
: 2Adam/m/conv_1/bias
: 2Adam/v/conv_1/bias
#:!	�. 2Adam/m/fc_0/kernel
#:!	�. 2Adam/v/fc_0/kernel
: 2Adam/m/fc_0/bias
: 2Adam/v/fc_0/bias
%:# 
2Adam/m/dense_3/kernel
%:# 
2Adam/v/dense_3/kernel
:
2Adam/m/dense_3/bias
:
2Adam/v/dense_3/bias
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
5__inference_random_color_affine_2_layer_call_fn_11635images"�
���
FullArgSpec!
args�
jimages

jtraining
varargs
 
varkw
 
defaults�
p

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
5__inference_random_color_affine_2_layer_call_fn_11640images"�
���
FullArgSpec!
args�
jimages

jtraining
varargs
 
varkw
 
defaults�
p

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
P__inference_random_color_affine_2_layer_call_and_return_conditional_losses_11679images"�
���
FullArgSpec!
args�
jimages

jtraining
varargs
 
varkw
 
defaults�
p

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
P__inference_random_color_affine_2_layer_call_and_return_conditional_losses_11683images"�
���
FullArgSpec!
args�
jimages

jtraining
varargs
 
varkw
 
defaults�
p

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
.__inference_gaussian_noise_layer_call_fn_11688inputs"�
���
FullArgSpec!
args�
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
.__inference_gaussian_noise_layer_call_fn_11693inputs"�
���
FullArgSpec!
args�
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
I__inference_gaussian_noise_layer_call_and_return_conditional_losses_11704inputs"�
���
FullArgSpec!
args�
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
I__inference_gaussian_noise_layer_call_and_return_conditional_losses_11708inputs"�
���
FullArgSpec!
args�
jinputs

jtraining
varargs
 
varkw
 
defaults�
p 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
v0
w1"
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
&__inference_conv_0_layer_call_fn_11717inputs"�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
A__inference_conv_0_layer_call_and_return_conditional_losses_11736inputs"�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
x0
y1"
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
&__inference_conv_1_layer_call_fn_11745inputs"�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
A__inference_conv_1_layer_call_and_return_conditional_losses_11764inputs"�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
'__inference_flatten_layer_call_fn_11769inputs"�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
B__inference_flatten_layer_call_and_return_conditional_losses_11775inputs"�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
z0
{1"
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
$__inference_fc_0_layer_call_fn_11784inputs"�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
?__inference_fc_0_layer_call_and_return_conditional_losses_11803inputs"�
���
FullArgSpec
args�

jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
__inference_loss_fn_0_11812"�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�B�
__inference_loss_fn_1_11821"�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�B�
__inference_loss_fn_2_11830"�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�B�
__inference_loss_fn_3_11839"�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�B�
__inference_loss_fn_4_11848"�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
�B�
__inference_loss_fn_5_11857"�
���
FullArgSpec
args� 
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *� 
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
 "
trackable_dict_wrapper�
 __inference__wrapped_model_10348w'()*+,%&8�5
.�+
)�&
input_5���������
'
� "1�.
,
dense_3!�
dense_3���������
�
A__inference_conv_0_layer_call_and_return_conditional_losses_11736s'(7�4
-�*
(�%
inputs���������
'
� "4�1
*�'
tensor_0���������# 
� �
&__inference_conv_0_layer_call_fn_11717h'(7�4
-�*
(�%
inputs���������
'
� ")�&
unknown���������# �
A__inference_conv_1_layer_call_and_return_conditional_losses_11764s)*7�4
-�*
(�%
inputs���������# 
� "4�1
*�'
tensor_0��������� 
� �
&__inference_conv_1_layer_call_fn_11745h)*7�4
-�*
(�%
inputs���������# 
� ")�&
unknown��������� �
B__inference_dense_3_layer_call_and_return_conditional_losses_11630c%&/�,
%�"
 �
inputs��������� 
� ",�)
"�
tensor_0���������

� �
'__inference_dense_3_layer_call_fn_11620X%&/�,
%�"
 �
inputs��������� 
� "!�
unknown���������
�
?__inference_fc_0_layer_call_and_return_conditional_losses_11803d+,0�-
&�#
!�
inputs����������.
� ",�)
"�
tensor_0��������� 
� �
$__inference_fc_0_layer_call_fn_11784Y+,0�-
&�#
!�
inputs����������.
� "!�
unknown��������� �
K__inference_finetuning_model_layer_call_and_return_conditional_losses_10884z'()*+,%&@�=
6�3
)�&
input_5���������
'
p

 
� ",�)
"�
tensor_0���������

� �
K__inference_finetuning_model_layer_call_and_return_conditional_losses_10931z'()*+,%&@�=
6�3
)�&
input_5���������
'
p 

 
� ",�)
"�
tensor_0���������

� �
K__inference_finetuning_model_layer_call_and_return_conditional_losses_11334y'()*+,%&?�<
5�2
(�%
inputs���������
'
p

 
� ",�)
"�
tensor_0���������

� �
K__inference_finetuning_model_layer_call_and_return_conditional_losses_11391y'()*+,%&?�<
5�2
(�%
inputs���������
'
p 

 
� ",�)
"�
tensor_0���������

� �
0__inference_finetuning_model_layer_call_fn_11000o'()*+,%&@�=
6�3
)�&
input_5���������
'
p

 
� "!�
unknown���������
�
0__inference_finetuning_model_layer_call_fn_11068o'()*+,%&@�=
6�3
)�&
input_5���������
'
p 

 
� "!�
unknown���������
�
0__inference_finetuning_model_layer_call_fn_11214n'()*+,%&?�<
5�2
(�%
inputs���������
'
p

 
� "!�
unknown���������
�
0__inference_finetuning_model_layer_call_fn_11235n'()*+,%&?�<
5�2
(�%
inputs���������
'
p 

 
� "!�
unknown���������
�
B__inference_flatten_layer_call_and_return_conditional_losses_11775h7�4
-�*
(�%
inputs��������� 
� "-�*
#� 
tensor_0����������.
� �
'__inference_flatten_layer_call_fn_11769]7�4
-�*
(�%
inputs��������� 
� ""�
unknown����������.�
I__inference_gaussian_noise_layer_call_and_return_conditional_losses_11704s;�8
1�.
(�%
inputs���������
'
p
� "4�1
*�'
tensor_0���������
'
� �
I__inference_gaussian_noise_layer_call_and_return_conditional_losses_11708s;�8
1�.
(�%
inputs���������
'
p 
� "4�1
*�'
tensor_0���������
'
� �
.__inference_gaussian_noise_layer_call_fn_11688h;�8
1�.
(�%
inputs���������
'
p
� ")�&
unknown���������
'�
.__inference_gaussian_noise_layer_call_fn_11693h;�8
1�.
(�%
inputs���������
'
p 
� ")�&
unknown���������
'C
__inference_loss_fn_0_11812$'�

� 
� "�
unknown C
__inference_loss_fn_1_11821$(�

� 
� "�
unknown C
__inference_loss_fn_2_11830$)�

� 
� "�
unknown C
__inference_loss_fn_3_11839$*�

� 
� "�
unknown C
__inference_loss_fn_4_11848$+�

� 
� "�
unknown C
__inference_loss_fn_5_11857$,�

� 
� "�
unknown �
P__inference_random_color_affine_2_layer_call_and_return_conditional_losses_11679s;�8
1�.
(�%
images���������
'
p
� "4�1
*�'
tensor_0���������
'
� �
P__inference_random_color_affine_2_layer_call_and_return_conditional_losses_11683s;�8
1�.
(�%
images���������
'
p 
� "4�1
*�'
tensor_0���������
'
� �
5__inference_random_color_affine_2_layer_call_fn_11635h;�8
1�.
(�%
images���������
'
p
� ")�&
unknown���������
'�
5__inference_random_color_affine_2_layer_call_fn_11640h;�8
1�.
(�%
images���������
'
p 
� ")�&
unknown���������
'�
G__inference_sequential_2_layer_call_and_return_conditional_losses_10570�'()*+,M�J
C�@
6�3
gaussian_noise_input���������
'
p

 
� ",�)
"�
tensor_0��������� 
� �
G__inference_sequential_2_layer_call_and_return_conditional_losses_10619�'()*+,M�J
C�@
6�3
gaussian_noise_input���������
'
p 

 
� ",�)
"�
tensor_0��������� 
� �
G__inference_sequential_2_layer_call_and_return_conditional_losses_11560w'()*+,?�<
5�2
(�%
inputs���������
'
p

 
� ",�)
"�
tensor_0��������� 
� �
G__inference_sequential_2_layer_call_and_return_conditional_losses_11611w'()*+,?�<
5�2
(�%
inputs���������
'
p 

 
� ",�)
"�
tensor_0��������� 
� �
,__inference_sequential_2_layer_call_fn_10682z'()*+,M�J
C�@
6�3
gaussian_noise_input���������
'
p

 
� "!�
unknown��������� �
,__inference_sequential_2_layer_call_fn_10744z'()*+,M�J
C�@
6�3
gaussian_noise_input���������
'
p 

 
� "!�
unknown��������� �
,__inference_sequential_2_layer_call_fn_11485l'()*+,?�<
5�2
(�%
inputs���������
'
p

 
� "!�
unknown��������� �
,__inference_sequential_2_layer_call_fn_11502l'()*+,?�<
5�2
(�%
inputs���������
'
p 

 
� "!�
unknown��������� �
G__inference_sequential_3_layer_call_and_return_conditional_losses_10394x@�=
6�3
)�&
input_6���������
'
p

 
� "4�1
*�'
tensor_0���������
'
� �
G__inference_sequential_3_layer_call_and_return_conditional_losses_10403x@�=
6�3
)�&
input_6���������
'
p 

 
� "4�1
*�'
tensor_0���������
'
� �
G__inference_sequential_3_layer_call_and_return_conditional_losses_11440w?�<
5�2
(�%
inputs���������
'
p

 
� "4�1
*�'
tensor_0���������
'
� �
G__inference_sequential_3_layer_call_and_return_conditional_losses_11444w?�<
5�2
(�%
inputs���������
'
p 

 
� "4�1
*�'
tensor_0���������
'
� �
,__inference_sequential_3_layer_call_fn_10414m@�=
6�3
)�&
input_6���������
'
p

 
� ")�&
unknown���������
'�
,__inference_sequential_3_layer_call_fn_10424m@�=
6�3
)�&
input_6���������
'
p 

 
� ")�&
unknown���������
'�
,__inference_sequential_3_layer_call_fn_11396l?�<
5�2
(�%
inputs���������
'
p

 
� ")�&
unknown���������
'�
,__inference_sequential_3_layer_call_fn_11401l?�<
5�2
(�%
inputs���������
'
p 

 
� ")�&
unknown���������
'�
#__inference_signature_wrapper_11169�'()*+,%&C�@
� 
9�6
4
input_5)�&
input_5���������
'"1�.
,
dense_3!�
dense_3���������
