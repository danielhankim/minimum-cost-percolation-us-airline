#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "binary_search_tree.h"

//code taken from https://www.programiz.com/dsa/binary-search-tree

// Create a node
struct node *newNode(int item) {
	struct node *temp = (struct node *)malloc(sizeof(struct node));
	temp->key = item;
	temp->left = temp->right = NULL;
	return temp;
}

// Inorder Traversal
void inorder(struct node *root) {
	if (root != NULL) {
		// Traverse left
		inorder(root->left);

		// Traverse root
		printf("%d -> ", root->key);

		// Traverse right
		inorder(root->right);
	}
}

// Insert a node
struct node *insert(struct node *node, int key) {
	// Return a new node if the tree is empty
	if (node == NULL) {return newNode(key);}

	// Traverse to the right place and insert the node
	if (key < node->key) {
		node->left = insert(node->left, key);
	} else {
		node->right = insert(node->right, key);
	}
	
	return node;
}

// Find the inorder successor
struct node *minValueNode(struct node *node) {
	struct node *current = node;

	// Find the leftmost leaf
	while (current && current->left != NULL) 
	current = current->left;

	return current;
}

// Deleting a node
struct node *deleteNode(struct node *root, int key) {
  // Return if the tree is empty
	if (root == NULL) {return root;}

	// Find the node to be deleted
	if (key < root->key) {
		root->left = deleteNode(root->left, key);
	} else if (key > root->key) {
		root->right = deleteNode(root->right, key);
	} else {
		// If the node is with only one child or no child
		if (root->left == NULL) {
			struct node *temp = root->right;
			free(root);
			return temp;
		} else if (root->right == NULL) {
			struct node *temp = root->left;
			free(root);
			return temp;
		}

		// If the node has two children
		struct node *temp = minValueNode(root->right);

		// Place the inorder successor in position of the node to be deleted
		root->key = temp->key;

		// Delete the inorder successor
		root->right = deleteNode(root->right, temp->key);
	}
	return root;
}



void free_tree(struct node *root){
	if (root != NULL) {
		free_tree(root->right);
		free_tree(root->left);
		free(root);
	}
	root = NULL;
}



// Inorder Traversal                                                                              
void tree_to_array (struct node *root, int *array) {
	if (root != NULL) {
		// Traverse left

		tree_to_array (root->left,array);


		// Traverse root                                                                              
		//printf("%d -> ", root->key);
		array[0]++;
		array[array[0]] = root->key;

		// Traverse right                                                                             
		tree_to_array (root->right,array);
	}
}
