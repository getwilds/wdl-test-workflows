version 1.0
## This is a test workflow that performs simple matrix 
## multiplication using a single GPU. 

#### WORKFLOW DEFINITION

workflow GpuMatrixMult {
  call Hello {
  }

  output {
    File stdout = Hello.response
  }

  parameter_meta {
    stdout: "printed results of the gpu script in question"
  }
}

#### TASK DEFINITIONS

task Hello {
  command <<<
    set -e
    echo 'hello there from the WILDS' 

    python3 <<CODE
    import tensorflow as tf
    # Creates a constant tensor from a tensor-like object.
    a = tf.constant([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], shape=[2, 3], name='a')
    b = tf.constant([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], shape=[3, 2], name='b')
    # Multiplies two tensors element-wise.
    c = tf.matmul(a, b)
    # Prints the result.
    print(c)
    # Check if GPU is available
    print("Num GPUs Available: ", len(tf.config.experimental.list_physical_devices('GPU')))
    CODE

    echo 'bye for now'
  >>>

  output {
    File response = stdout()
  }

  runtime {
    docker: 'tensorflow/tensorflow:latest-gpu'
    # modules: "TensorFlow/2.11.0-foss-2022a-CUDA-11.7.0"
    gpus: '1'
  }

  parameter_meta {
    response: "printed results of the gpu script in question"
  }
}

