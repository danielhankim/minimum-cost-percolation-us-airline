struct node *newNode(int item);
void inorder(struct node *root);
struct node *insert(struct node *node, int key);
struct node *minValueNode(struct node *node);
struct node *deleteNode(struct node *root, int key);
void free_tree(struct node *root);
void tree_to_array(struct node *root, int *array);
