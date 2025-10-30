import torch


def small_to_unknown(sig_weights_batch, cutoff=0.01):
    """Send all values < cutoff to 0, adding their sum to a new element

    Args:
        sig_weights_batch (torch.tensor): Batch of signature weights to normalie
        cutoff (float): Value from which to send to zero
    """
    mask = (sig_weights_batch > cutoff).to(torch.int)
    sig_weights_batch = mask*sig_weights_batch
    result = torch.zeros((sig_weights_batch.shape[0], sig_weights_batch.shape[1] + 1)).to(sig_weights_batch)
    result[:, :-1] = sig_weights_batch
    result[:, -1] = (1 - torch.sum(sig_weights_batch, axis=1))
    return result


if __name__ == "__main__":
    example = torch.tensor([
        [0.1, 0.7, 0.2],
        [0.05, 0.1, 0.85],
        [0.5, 0.5, 0],
        [0.5, 0.5, 0],
    ])

    cutoff = 0.15

    example_with_unknown = small_to_unknown(example, cutoff)

    print(example)
    print(example_with_unknown)
